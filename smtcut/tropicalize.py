#! /usr/bin/env python3

# $+HEADER$
#
# Copyright 2019-2021 Christoph Lueders
#
# This file is part of the SMTcut project: <http://wrogn.com/smtcut>
#
# SMTcut is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SMTcut is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SMTcut.  If not, see <http://www.gnu.org/licenses/>.
#
# $-HEADER$

# is imported by: smtcut

# pylance local imports issue:
# https://github.com/microsoft/pylance-release/issues/68

import sys, os, re
from sympy import sympify, Rational, sstr, gcd, Symbol, lcm_list, Mul
from math import log
from operator import mod
from gmpy2 import mpq
from pprint import pprint

from .util import mytime
from .phkeep import Bag, PhKeep, mpq_str


def nstr(f):
    return sstr(f, order="none")


class CollectProd:
    """
    Collect a product of parameters (and maybe a literal factor).
    """
    def __init__(self, log):
        self.prod = 1
        self.log = log
    def __call__(self, x=None, pow=1):
        if x is not None:
            self.prod *= x ** pow
    def get(self):
        return self.log(self.prod)

class CollectSum:
    """
    Collect a sum of rounded logarithms of parameters (and maybe a literal factor).
    """
    def __init__(self, log):
        self.sum = 0
        self.log = log
    def __call__(self, x=None, pow=1):
        if x is not None:
            self.sum += pow * self.log(x)
    def get(self):
        return self.sum


def make_mpq(i):
    """
    Convert to gmpy2.mpq.  Accepts int, Rational

    >>> make_mpq(-2)
    mpq(-2,1)
    >>> make_mpq(Rational(2,3))
    mpq(2,3)
    >>> make_mpq(0.25)
    mpq(1,4)
    """
    try:
        # int, float, string
        return mpq(i)
    except TypeError:
        # Rational
        return mpq(i.numerator(), i.denominator())


def _gcd_tuple(d):
    g = make_mpq(gcd(d))
    if g == 0:
        return tuple(d)
    if g < 0:
        # must not change sign
        g = -g
    return tuple(i / g for i in d)


split_pattern = re.compile(r"\A(?:([a-zA-Z]+)_?([0-9]+)|([a-zA-Z]+))\Z")
inf = float("inf")

def split_name_num(s):
    """
    Split a variable name into "base name" and index.
    Index will be -inf, if no index was given.

    >>> split_name_num("x1")
    ('x', 1)
    >>> split_name_num("x_1")
    ('x', 1)
    >>> split_name_num("x0")
    ('x', 0)
    >>> split_name_num("x10")
    ('x', 10)
    >>> split_name_num(sympify("x10"))
    ('x', 10)
    >>> split_name_num("x")
    ('x', -inf)
    >>> split_name_num("$") is None
    True
    >>> split_name_num("x") < split_name_num("x0")
    True
    """
    s = str(s)
    m = split_pattern.match(s)
    if not m:
        return None
    if m.group(3):
        return m.group(3), -inf
    return m.group(1), int(m.group(2))


def logep(gamma, eps, p=Rational(1)):
    """
    Compute a c from gamma's absolute value.

    >>> logep(2, Rational(1,2), 1)
    -1
    >>> logep(Rational("1.3"), Rational(1,2), 10)
    -2/5
    >>> logep(Rational("1.3"), Rational(1,2), 100)
    -19/50
    """
    assert p == int(p) and p > 0
    p = Rational(p)
    assert gamma != 0
    return round(p * log(gamma * (-1 if gamma < 0 else 1), eps)) / p


class Tropicalize:
    def __init__(self, cfg, *, evaluate_literals=False, round_once=False, substitute_first=False, param={}, vars=[], slows=[]):
        """
        logep: log function with round and scale to base epsilon in (0,1).
        param: dict with parameter names (str) and values (number).
        vars: a list of variables (str) that specifies their order.
        """
        self.cfg = cfg
        self.eps = cfg.eps
        self.p = cfg.p
        self.prt = cfg.prt
        self.verbose = cfg.verbose
        self.evaluate_literals = evaluate_literals
        self.round_once = round_once
        self.substitute_first = substitute_first
        self.param = param
        self.eqns = []
        self.vars = vars
        self.slow = slows
        self.size_vector_field = 0

    def _load_parameters(self, seq, fn):
        """
        >>> t = Tropicalize()
        >>> seq = ["k1 = 12", "k2=1.5", " k3 =  -12/13", "k4 = -1e-75"]
        >>> t._load_parameters(seq, "-")
        {'k1': 12, 'k2': 3/2, 'k3': -12/13, 'k4': -1/1000000000000000000000000000000000000000000000000000000000000000000000000000}
        """
        d = {}
        line = 0

        for ln in seq:
            line += 1
            ln = ln.strip()
            # ignore empty or comment lines
            if ln[:1] in ("", ";", "#"):
                continue
            a = ln.split("=")
            if len(a) != 2:
                self.prt(f"{fn} ({line}): parameter lines must be: VARIABLE = VALUE\n", screen=self.verbose >= 0)
                continue
            d[a[0].strip()] = Rational(a[1].strip())

        return d

    def load_parameters(self, fn):
        """
        Read file 'fn' as key=value pairs and save a dict with str:Rational pairs as self.param.
        """
        self.prt(f"Loading parameters from {fn} ...", screen=self.verbose >= 1)
        with open(fn) as seq:
            self.param = self._load_parameters(seq, fn)

    def _load_equations(self, seq):
        """
        >>> t = Tropicalize()
        >>> seq = ["k2*x1 + 2*x2 - x1*x2", "x2**k1 - log(3)*x3", "x + z"]
        >>> print(nstr(t._load_equations(seq)))
        ([k2*x1 + 2*x2 - x1*x2, x2**k1 - log(3)*x3, x + z], [k1, k2, x, x1, x2, x3, z])
        """
        eqns = []
        line = 0
        vars = set()

        for ln in seq:
            ln = ln.strip()
            # ignore empty or comment lines
            if ln[:1] in ("", ";", "#"):
                continue
            f = sympify(ln, evaluate=False)
            eqns.append(f)
            vars.update(str(i) for i in f.free_symbols)

        # build list of variables, filter out parameter names and sort properly
        vars = sorted(vars - set(self.param.keys()), key=lambda x: split_name_num(x))
        return eqns, vars

    def load_equations(self, fns):
        """
        Read file 'fns' as r.h.s. of 0 = f and
        return a list of sympy expressions and a sorted list of variables.
        Also accepts list of filenames, which are then read one after the other to form one set of input equations.
        """
        l = []
        if not any(isinstance(fns, i) for i in (list, tuple)):
            fns = [fns]
        for fn in fns:
            self.prt(f"Loading {fn} ...", screen=self.verbose >= 1)
            with open(fn) as seq:
                l.extend(seq.readlines())
                if fn == fns[0]:
                    self.size_vector_field = len(l) #todo: it remains a problem if comment lines
        self.eqns, self.vars = self._load_equations(l)

    def _load_slows(self, seq, f):
        """
        >>> t = Tropicalize()
        >>> seq = ["x2", "x5"]
        >>> t._load_slows(seq, "-")
        ['x2','x5']
        """
        v = []
        line = 0
        for ln in seq:
            ln = ln.strip()
            # ignore empty or comment lines
            if ln[:1] in ("", ";", "#"):
                continue
            # todo : add verification
            v.append(ln)
        return v

    def load_slows(self, fn):
        self.prt(f"Loading {fn} ...", screen=self.verbose >= 1)
        with open(fn) as seq:
            self.slow = self._load_slows(seq, fn)

    def T_trop_term(self, f):
        """
        >>> t = Tropicalize(2, evaluate_literals=True, vars=["x", "y", "z"])
        >>> print(t.T_trop_term(sympify("x**2+2*y-7*z**5")))
        ([(0, 2, 0, 0), (1, 0, 1, 0)], [(3, 0, 0, 5)])

        >>> t = Tropicalize(2, param={"k1":2, "k2":7}, vars=["x", "y", "z"])
        >>> print(t.T_trop_term(sympify("x**2+k1*y-k2*z**5")))
        ([(0, 2, 0, 0), (1, 0, 1, 0)], [(3, 0, 0, 5)])
        
        >>> t = Tropicalize(2, param={"k1":2, "k2":2}, vars=["x", "y"])
        >>> print(t.T_trop_term(sympify("x**2+k1*y-k2*y")))
        ([(0, 2, 0), (1, 0, 1)], [(1, 0, 1)])
        
        >>> t = Tropicalize(2, evaluate_literals=True, substitute_first=True, param={"k1":2, "k2":2}, vars=["x", "y"])
        >>> print(t.T_trop_term(sympify("x**2+k1*y-k2*y")))
        ([(0, 2, 0)], [])
        """
        pp, np = self.trop_term(f)
        return mpq_str((sorted(pp), sorted(np)))

    def trop_term(self, f, eqcnt):
        """
        Tropicalize equation f = 0, where f is a sympy expression.
        Returns (positive points, negative points), where the coordinates are of type 'mpq'.
        """
        dim = len(self.vars)
        assert len(set(self.vars)) == dim                        # no double use
        assert not (set(self.vars) & set(self.param.keys()))          # must be disjunct
        assert not self.substitute_first or self.evaluate_literals    # substitute_first requires not evaluate_literals
        assert set(str(i) for i in f.free_symbols) <= set(self.vars) | set(self.param.keys())  # no unknown symbols

        # evaluate to get rid of log(2) and the like
        if self.verbose >= 3:
            self.prt(f"\nInput:    {f}")
        if self.substitute_first:
            f = f.subs(self.param, evaluate=False)
        f = f.evalf().expand()
        if self.verbose >= 3:
            self.prt(f"Expanded: {f}")
        if f == 0:
            self.prt(end=f"\n*** Equation {eqcnt} is zero, ignoring. ", screen=self.verbose >= 2)
            return None
        tt, v = f.as_terms()

        # check for denominators
        denoms = []
        for _,((re,im),monom,ncpart) in tt:
            for i,m in enumerate(monom):
                # only consider denominators that contain variables
                if m < 0 and type(v[i]) != Symbol and set(str(i) for i in v[i].free_symbols) - set(self.param.keys()):
                    # denominator found
                    lfact = []
                    for fact in v[i].factor().as_ordered_factors():
                        if set(str(i) for i in fact.free_symbols) - set(self.param.keys()):
                            # this factor contains variables
                            lfact.append(fact)
                    if lfact:
                        denoms.append(Mul(*lfact) ** -m)

        if denoms:
            # multiply equation with lcm of all denoms
            denom = lcm_list(denoms).factor()
            self.prt(end=f"\n*** Multiplying equation {eqcnt} with {denom}. ", screen=self.verbose >= 1)
            f = (f * denom).cancel().expand()
            if self.verbose >= 3:
                self.prt(f"Canceled: {f}")
            tt, v = f.as_terms()

        assert all(type(i) == Symbol for i in v)            # all 'variables' in 'v' must be symbols
        v = [str(i) for i in v]                             # convert to str
        assert set(v) <= set(self.vars) | set(self.param.keys())  # must be listed in self.vars

        # build translation list between self.vars and v (mm), and
        # list of parameter values (pv)
        mm = [-1] * len(v)
        pv = [0] * len(v)
        for i,s in enumerate(v):
            try:
                mm[i] = self.vars.index(str(s)) + 1
            except ValueError:
                pv[i] = self.param[str(s)]

        # set of positive and negative points
        pp, np = set(), set()
        mylog = lambda x: logep(x, self.eps, self.p)

        # cycle through all monomials
        for term,((re,im),monom,ncpart) in tt:
            # each tuple is (term, ((real, imag), monom, ncpart))
            assert re != 0
            assert im == 0
            assert not ncpart
            # print(f"term={_}")

            sign = 1 if re > 0 else -1
            # either add the rounded logs, or multiply the values and do the rounded log once at the end
            collect = (CollectProd if self.round_once else CollectSum)(mylog)
            p = [0] * (1 + dim)
            if self.evaluate_literals:
                # in tropicalization, usually literals are ignored: i.e., Trop(2*x) = Trop(x)
                collect(sign * re)
            #print(monom, v)
            for i,m in enumerate(monom):
                idx = mm[i]
                if idx == -1:
                    # parameter
                    assert not self.substitute_first        # can not happen otherwise
                    # print(m, pv[i])
                    if pv[i] == 0:
                        if m != 0:
                            self.prt(end=f"\n*** Monomial {term} in equation {eqcnt} is zero, since {v[i]}=0. ", screen=self.verbose >= 2)
                            sign = 0
                            break
                    else:
                        collect(pv[i], m)
                else:
                    # variable
                    p[idx] = m
            # save the value of the rounded log collection
            p[0] = collect.get()

            # add point to the right set
            if sign != 0:
                (pp if sign == 1 else np).add((tuple(make_mpq(i) for i in p),(term,eqcnt)))

        if not pp and not np:
            self.prt(end=f"\n*** Equation {eqcnt} has no monomials, ignoring. ", screen=self.verbose >= 1)
            return None

        if not pp or not np:
            self.prt(end=f"\n*** Equation {eqcnt} has unbalanced monomials, ignoring. ", screen=self.verbose >= 1)
            return None

        return pp, np

    def trop_term_partial(self, f1, g1, eqcnt, eqcntb, scnt):

        dim = len(self.vars)
        assert len(set(self.vars)) == dim                        # no double use
        assert not (set(self.vars) & set(self.param.keys()))          # must be disjunct
        assert not self.substitute_first or self.evaluate_literals    # substitute_first requires not evaluate_literals
        assert set(str(i) for i in f1.free_symbols) <= set(self.vars) | set(self.param.keys())  # no unknown symbols
        assert set(str(i) for i in g1.free_symbols) <= set(self.vars) | set(self.param.keys())

        #symvar = sympify(self.vars)
        varslow = 'x'+ str(eqcnt)
        varfast = 'x' + str(eqcntb)
        symvarslow = sympify(varslow)
        symvarfast = sympify(varfast)
        f = f1/symvarslow
        g = g1/symvarfast

        # evaluate to get rid of log(2) and the like
        if self.verbose >= 3:
            self.prt(f"\nInput:    {f}")
        if self.substitute_first:
            f = f.subs(self.param, evaluate=False)
        f = f.evalf().expand()
        if self.verbose >= 3:
            self.prt(f"Expanded: {f}")
        if f == 0:
            self.prt(end=f"\n*** Equation {eqcnt} is zero, ignoring. ", screen=self.verbose >= 2)
            return None
        ttf, vf = f.as_terms()

        # evaluate to get rid of log(2) and the like
        if self.verbose >= 3:
            self.prt(f"\nInput:    {g}")
        if self.substitute_first:
            g = g.subs(self.param, evaluate=False)
        g = g.evalf().expand()
        if self.verbose >= 3:
            self.prt(f"Expanded: {g}")
        if g == 0:
            self.prt(end=f"\n*** Equation {eqcntb} is zero, ignoring. ", screen=self.verbose >= 2)
            return None
        ttg, vg = g.as_terms()

        # check for denominators
        denomsf = []
        for _,((re,im),monom,ncpart) in ttf:
            for i,m in enumerate(monom):
                # only consider denominators that contain variables
                if m < 0 and type(vf[i]) != Symbol and set(str(i) for i in vf[i].free_symbols) - set(self.param.keys()):
                    # denominator found
                    lfact = []
                    for fact in vf[i].factor().as_ordered_factors():
                        if set(str(i) for i in fact.free_symbols) - set(self.param.keys()):
                            # this factor contains variables
                            lfact.append(fact)
                    if lfact:
                        denomsf.append(Mul(*lfact) ** -m)

        if denomsf:
            # multiply equation with lcm of all denoms
            denomf = lcm_list(denomsf).factor()
            self.prt(end=f"\n*** Multiplying equation {eqcnt} with {denomf}. ", screen=self.verbose >= 1)
            f = (f * denomf).cancel().expand()
            if self.verbose >= 3:
                self.prt(f"Canceled: {f}")
            ttf, vf = f.as_terms()

        # check for denominators
        denomsg = []
        for _,((re,im),monom,ncpart) in ttg:
            for i,m in enumerate(monom):
                # only consider denominators that contain variables
                if m < 0 and type(vg[i]) != Symbol and set(str(i) for i in vg[i].free_symbols) - set(self.param.keys()):
                    # denominator found
                    lfact = []
                    for fact in vg[i].factor().as_ordered_factors():
                        if set(str(i) for i in fact.free_symbols) - set(self.param.keys()):
                            # this factor contains variables
                            lfact.append(fact)
                    if lfact:
                        denomsg.append(Mul(*lfact) ** -m)

        if denomsg:
            # multiply equation with lcm of all denoms
            denomg = lcm_list(denomsg).factor()
            self.prt(end=f"\n*** Multiplying equation {eqcnt} with {denomg}. ", screen=self.verbose >= 1)
            g = (g * denomg).cancel().expand()
            if self.verbose >= 3:
                self.prt(f"Canceled: {g}")
            ttg, vg = g.as_terms()

        assert all(type(i) == Symbol for i in vf)            # all 'variables' in 'v' must be symbols
        vf = [str(i) for i in vf]                             # convert to str
        assert set(vf) <= set(self.vars) | set(self.param.keys())  # must be listed in self.vars
        
        assert all(type(i) == Symbol for i in vg)            # all 'variables' in 'v' must be symbols
        vg = [str(i) for i in vg]                             # convert to str
        assert set(vg) <= set(self.vars) | set(self.param.keys())  # must be listed in self.vars

        # build translation list between self.vars and v (mm), and
        # list of parameter values (pv)
        mmf = [-1] * len(vf)
        pvf = [0] * len(vf)
        for i,s in enumerate(vf):
            try:
                mmf[i] = self.vars.index(str(s)) + 1
            except ValueError:
                pvf[i] = self.param[str(s)]

        # set of slow points
        sp = set()
        mylogf = lambda x: logep(x, self.eps, self.p)

         # cycle through all monomials
        for term,((re,im),monom,ncpart) in ttf:
            # each tuple is (term, ((real, imag), monom, ncpart))
            assert re != 0
            assert im == 0
            assert not ncpart
            # print(f"term={_}")

            sign = 1 if re > 0 else -1
            # either add the rounded logs, or multiply the values and do the rounded log once at the end
            collect = (CollectProd if self.round_once else CollectSum)(mylogf)
            p = [0] * (1 + dim)
            if self.evaluate_literals:
                # in tropicalization, usually literals are ignored: i.e., Trop(2*x) = Trop(x)
                collect(sign * re)
            #print(monom, v)
            for i,m in enumerate(monom):
                idx = mmf[i]
                if idx == -1:
                    # parameter
                    assert not self.substitute_first        # can not happen otherwise
                    # print(m, pvf[i])
                    if pvf[i] == 0:
                        if m != 0:
                            self.prt(end=f"\n*** Monomial {term} in equation {eqcnt} is zero, since {vf[i]}=0. ", screen=self.verbose >= 2)
                            sign = 0
                            break
                    else:
                        collect(pvf[i], m)
                else:
                    # variable
                    p[idx] = m
            # save the value of the rounded log collection
            p[0] = collect.get()

            sp.add((tuple(make_mpq(i) for i in p),(term,eqcnt)))

        if not sp:
            self.prt(end=f"\n*** Equation {eqcnt} has no monomials, ignoring. ", screen=self.verbose >= 1)
            return None

        # build translation list between self.vars and v (mm), and
        # list of parameter values (pv)
        mmg = [-1] * len(vg)
        pvg = [0] * len(vg)
        for i,s in enumerate(vg):
            try:
                mmg[i] = self.vars.index(str(s)) + 1
            except ValueError:
                pvg[i] = self.param[str(s)]

        # set of fast points
        fp = set()
        mylogg = lambda x: logep(x, self.eps, self.p)

         # cycle through all monomials
        for term,((re,im),monom,ncpart) in ttg:
            # each tuple is (term, ((real, imag), monom, ncpart))
            assert re != 0
            assert im == 0
            assert not ncpart
            # print(f"term={_}")

            sign = 1 if re > 0 else -1
            # either add the rounded logs, or multiply the values and do the rounded log once at the end
            collect = (CollectProd if self.round_once else CollectSum)(mylogg)
            p = [0] * (1 + dim)
            if self.evaluate_literals:
                # in tropicalization, usually literals are ignored: i.e., Trop(2*x) = Trop(x)
                collect(sign * re)
            #print(monom, v)
            for i,m in enumerate(monom):
                idx = mmg[i]
                if idx == -1:
                    # parameter
                    assert not self.substitute_first        # can not happen otherwise
                    # print(m, pvg[i])
                    if pvg[i] == 0:
                        if m != 0:
                            self.prt(end=f"\n*** Monomial {term} in equation {eqcntb} is zero, since {vg[i]}=0. ", screen=self.verbose >= 2)
                            sign = 0
                            break
                    else:
                        collect(pvg[i], m)
                else:
                    # variable
                    p[idx] = m
            # save the value of the rounded log collection
            p[0] = collect.get()

            fp.add((tuple(make_mpq(i) for i in p),(term,eqcntb)))

        if not fp:
            self.prt(end=f"\n*** Equation {eqcnt} has no monomials, ignoring. ", screen=self.verbose >= 1)
            return None

        return sp, fp

    def build_polyhedra(self, pp, np):
        b = Bag()
        from itertools import product

        # prt(end=f"{len(pp),len(np)}", flush=True)
        for v1, v2 in product(pp, np):
            # v1 = v2, hence v1 - v2 = 0.  This defines a hyperplane.
            d = _gcd_tuple([a - b for a, b in zip(v1[0], v2[0])])
            # Furthermore, v1 = v2 <= vi (where vi \in pp \cup np),
            # hence vi - v1 >= 0.  This defines the half-spaces.
            ie = []
            for vi in set(pp) | set(np):
                if vi != v1 and vi != v2:
                    ie.append(_gcd_tuple([a - b for a, b in zip(vi[0], v1[0])]))
            # create the polyhedron
            p = PhKeep(eqns=[d], ieqs=ie, terms=[(v1[1],v2[1])])  
            if not p.is_empty(self.cfg):
                b.insert_include(self.cfg, p.minimize(self.cfg))

        return b

    def build_polyhedra_partial(self, sp, fp):
        b = Bag()
        from itertools import product

        #sp represent slow_points and fp fast-points
        for v1, v2 in product(sp, fp):
            ie = []
            #v2 <= v1, hence v1 - v2 >= 0.
            ie.append(_gcd_tuple([a - b for a, b in zip(v1[0], v2[0])]))
            for v3 in set(sp):
                if v3 != v1:
                    ie.append(_gcd_tuple([b - a for a, b in zip(v1[0], v3[0])]))
            for v4 in set(fp):
                if v4 != v2:
                    ie.append(_gcd_tuple([b - a for a, b in zip(v2[0], v4[0])]))
            p = PhKeep(eqns=[], ieqs=ie, terms=[(v1[1],v2[1])])
            if not p.is_empty(self.cfg):
                b.insert_include(self.cfg, p.minimize(self.cfg))

        return b

    def tropicalize(self, param_fn, eqn_fns, slow_fn, ret=None):
        """
        From 'dir' load equations and tropicalize them with epsilon='self.eps' and p='self.p'.
        """
        self.load_parameters(param_fn)
        # print(self.param)
        self.load_equations(eqn_fns)
        # pprint(self.eqns)
        self.load_slows(slow_fn)
        self.slow.sort()
        #todo: add verification about slow

        self.prt(end=f"Tropicalizing (eps={self.eps}, p={self.p}) ... ", screen=self.verbose >= 1)
        t = mytime()
        bb = []
        if len(self.slow)==0:
            for cnt, f in enumerate(self.eqns, 1):
                self.prt(end=f"[{cnt}]", flush=True, screen=self.verbose >= 1)
                r = self.trop_term(f, cnt)
                if r is not None:
                    b = self.build_polyhedra(*r)
                    bb.append(b)
        else:
        #the first loop is used to get equations linked to fast species.
            for cnt, f in enumerate(self.eqns, 1):
                fast = True
                for c, i in enumerate(self.slow):
                    if cnt == split_name_num(self.slow[c])[1]:
                        fast = False
                        #break
                if fast:#if f is fast
                    self.prt(end=f"[{cnt}]", flush=True, screen=self.verbose >= 1)
                    r0 = self.trop_term(f, cnt)#equilibration
                    if r0 is not None:
                        b = self.build_polyhedra(*r0)
                        bb.append(b)
                    if cnt<=self.size_vector_field:
                        cnt_slow = 0
                        for cntb, g in enumerate(self.eqns, 1):
                            if cnt_slow<len(self.slow):
                                nn = split_name_num(self.slow[cnt_slow])
                                if cntb == nn[1]:#if g is slow
                                    self.prt(end=f"[{cntb}]", flush=True, screen=self.verbose >= 1)
                                    r1 = self.trop_term_partial(g, f, cntb, cnt, cnt_slow)#slow/fast splitting
                                    if r1 is not None:
                                        b = self.build_polyhedra_partial(*r1)
                                        bb.append(b)
                                    cnt_slow+=1
        t = mytime() - t
        self.prt(f"\nTropicalization {t:.3f} sec\n", flushfile=True, screen=self.verbose >= 1)

        if ret:
            ret.tropicalization_time = t
                    
        return bb


def tester():
    from .prt import prt
    prt("test: tropicalize")
    import doctest, sys
    this_mod = sys.modules[__name__]    
    doctest.testmod(this_mod, verbose=False)
