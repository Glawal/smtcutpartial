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

# is imported by: preproc, main, smtcut, tropicalize

import re
from operator import mul

from .smt import *


_mpq_str_pat = re.compile(r"mpq\((-?[0-9]+),([0-9]+)\)")
def _mpq_repl(m):
    return m.group(1) if m.group(2) == "1" else f"{m.group(1)}/{m.group(2)}"

def mpq_str(i):
    """
    Print nested structures with "mpq()" around rationals.

    >>> print(mpq_str(mpq(2,3)))
    2/3
    >>> print(mpq_str([mpq(2,-3), mpq(4,-2)]))
    [-2/3, -2]
    >>> print(mpq_str((mpq(2,3), mpq(-2,1))))
    (2/3, -2)
    """
    s = str(i)
    return _mpq_str_pat.sub(_mpq_repl, s)


class Bag:
    """
    A bag is a list of polyhedra.
    """
    def __init__(self, l=None, name=None):
        self.l = l if l is not None else []
        self.name = "" if name is None else name

    def __len__(self):
        return len(self.l)

    def __getitem__(self, length):
        return self.l[length]

    def __setitem__(self, key, value):
        self.l[key] = value

    def append(self, n):
        self.l.append(n)

    def extend(self, n):
        self.l.extend(n)

    def __str__(self):
        return f'[{", ".join([str(i) for i in self.l])}]'

    def insert_larger(self, cfg, n):
        """
        Insert polyhedron n into list self.l.  We assume that n is not included
        in any of self.l's members.  Check if n includes any of self.l's members.
        """
        self._remove_included(cfg, n, [n])

    def remove_included(self, cfg, n):
        """
        Remove all polyhedra from self.l that are included in polyhedron n.
        """
        self._remove_included(cfg, n, [])

    def _remove_included(self, cfg, n, l2):
        """
        Remove all polyhedra from self.l that are included in polyhedron n.
        """
        for p in self.l:
            # is p contained in n?  if so, don't copy it to new list.
            if not n.contains(cfg, p):
                # append p to new list
                l2.append(p)
        self.l = l2

    def insert_include(self, cfg, n):
        """
        Insert polyhedron n into list self.l, but check that n is not included
        in any of self.l's members and that n doesn't include any of self.l's members.
        If so, use the larger one.
        """
        # we're assuming that it is more likely that the new polyhedron is already included in
        # some polyhedron already in the list, so we check for that first.
        for p in self.l:
            # is n contained in p, i.e. something we already know?
            if p.contains(cfg, n):
                return                                      # no need to continue
        # n is not included in any of l's polyhedra, so it has to be added to the list.
        # see if n includes any of the already existing ones.
        l2 = [n]
        for p in self.l:
            # is p contained in n?  if so, don't copy it to new list.
            if not n.contains(cfg, p):
                l2.append(p)                                # append p to new list
        self.l = l2


class PhKeep:
    def __init__(self, eqns=None, ieqs=None, phkeep=None, terms=None):
        self.eqns = set()
        self.ieqs = set()
        self.terms = []
        if eqns:
            self.eqns.update(eqns)
        if ieqs:
            self.ieqs.update(ieqs)
        if terms:
            self.terms += terms
        if phkeep is not None:
            self.eqns.update(phkeep.eqns)
            self.ieqs.update(phkeep.ieqs)
            self.terms += phkeep.terms
        #self.point = None
        self.f = phkeep.f if phkeep and not ieqs and not eqns else None
        

    def __str__(self):
        # remove spaces in tuples of (in)equalities; sort sets for printing
        return f"<eqns=[{', '.join([mpq_str(i).replace(' ','') for i in sorted(self.eqns)])}] ieqs=[{', '.join([mpq_str(i).replace(' ','') for i in sorted(self.ieqs)])}]>"

    def add_constraints(self, *, eqns=None, ieqs=None, phkeep=None, terms=None):
        if eqns:
            self.eqns.update(eqns)
        if ieqs:
            self.ieqs.update(ieqs)
        if terms:
            self.terms += terms
        if phkeep is not None:
            self.eqns.update(phkeep.eqns)
            self.ieqs.update(phkeep.ieqs)
            self.terms += phkeep.terms
        self.f = None

    def is_universe(self):
        return self.eqns is None and self.ieqs is None

    def _make_f(self, cfg):
        """
        Create a formula for the polyhedron and save it for later.
        If such formula was already saved, use it.
        """
        if self.f is None:
            self.f = cfg.mangler.smt_ph(self)
        return self.f

    def test_empty(self, cfg):
        with cfg.aux_solver as saux:
            r = saux.is_sat(self._make_f(cfg))
        return not r

    def test_intersection_empty(self, cfg, other):
        with cfg.aux_solver as saux:
            r = saux.is_sat(And(self._make_f(cfg), other._make_f(cfg)))
        return not r

    def intersect(self, other):
        """
        Just add the two sets of constraints.  Doesn't perform any minimization.
        """
        phnew = PhKeep(phkeep=self)
        phnew.add_constraints(phkeep=other)
        return phnew

    def space_dim(self):
        """
        Dimension of the ambient space.
        """
        return len(next(iter(self.eqns|self.ieqs))) - 1

    def contains_point(self, v):
        """
        Does self contain point v?

        >>> PhKeep(eqns=[(-1,1,-1)], ieqs=[(-2,1,0)]).contains_point((3,2))
        True
        >>> PhKeep(eqns=[(-1,1,-1)], ieqs=[(-2,1,0)]).contains_point((3,3))
        False
        """
        # obvious version, fastest
        return (all(sum(map(mul, v, eq[1:])) == -eq[0] for eq in self.eqns) and
            all(sum(map(mul, v, ie[1:])) >= -ie[0] for ie in self.ieqs))

    def contains_smt(self, cfg, other):
        """
        Does self contain polyhedron other?  Full SMT check.
        """
        f = And(other._make_f(cfg), Not(self._make_f(cfg)))
        with cfg.aux_solver as saux:
            r = saux.is_sat(f)
        return not r

    def contains(self, cfg, other):
        """
        Does self contain polyhedron other?
        Uses other.point to weed out most polyhedra that don't include the point.
        """
        try:
            if not self.contains_point(other.point):
                return False
        except AttributeError:
            # 'other' has no attribute 'point'
            pass
        return self.contains_smt(cfg, other)

    def minimize(self, cfg, *, min_eqns=None, i2e=False):
        """
        Go through all constraints and use SMT to test if they are required.
        """
        if not self.eqns and not self.ieqs:
            return

        nieqs = self.ieqs 
        self.ieqs = set() 

        if min_eqns is None:
            neqns = self.eqns 
            self.eqns = set() 
        else:
            # min_eqns must be a subset of self.eqns
            assert min_eqns <= self.eqns
            
            # filter out already minimized equalities (!) only
            neqns = self.eqns - min_eqns
            self.eqns = set(min_eqns)       # must copy

        self.f = None

        # sorting the (in)equalities by number of non-zero coefficients
        # yields some 6% speed-up. (BM103 -e11 -msat)

        with cfg.aux_solver as saux:
            if min_eqns is not None:
                saux.add_assertion(cfg.mangler.smt_eqns(min_eqns))

            # walk through all assertions
            for q in sorted(neqns, key=lambda x: sum(i != 0 for i in x[1:])):
            # for q in neqns:
                c = cfg.mangler.smt_eqn(q)
                saux.push()
                r = saux.is_sat(Not(c))
                saux.pop()
                if r:
                    # g is becoming smaller, i.e., q was needed
                    saux.add_assertion(c)
                    self.eqns.add(q)

            for q in sorted(nieqs, key=lambda x: sum(i != 0 for i in x[1:])):
            # for q in nieqs:
                c = cfg.mangler.smt_ieq(q)
                saux.push()
                r = saux.is_sat(Not(c))
                saux.pop()
                if r:
                    # g is becoming smaller, i.e., q was needed
                    saux.add_assertion(c)
                    self.ieqs.add(q)

            # THIS DOESN'T WORK PROPERLY YET!!!
            if i2e:
                # check if half-spaces contribute to hyperplanes
                nieqs = self.ieqs
                self.ieqs = set()
                for q in nieqs:
                    # treat ieqn as eqn and see if it changes anything
                    c = cfg.mangler.smt_eqn(q)
                    saux.push()
                    r = saux.is_sat(Not(c))
                    saux.pop()
                    if not r:
                        # g doesn't become smaller, i.e. adding the ieqn as eqn does not change the definition; keep q as equation
                        self.eqns.add(q)
                    else:
                        self.ieqs.add(q)

                # re-check which ieqs are needed
                if len(nieqs) != len(self.ieqs):
                    nieqs = self.ieqs
                    self.ieqs = set()
                    for q in nieqs:
                        c = cfg.mangler.smt_ieq(q)
                        saux.push()
                        r = saux.is_sat(Not(c))
                        saux.pop()
                        if r:
                            # g is becoming smaller, i.e., q was needed
                            saux.add_assertion(c)
                            self.ieqs.add(q)

        return self

    def is_empty(self, cfg):
        """
        Does any point satisfy the constraints?
        >>> PhKeep(eqns=[], ieqs=[]).is_empty()
        False
        >>> PhKeep(eqns=[], ieqs=[(0, 1, 1)]).is_empty()
        False
        >>> PhKeep(eqns=[], ieqs=[(-1, 1), (0, -1)]).is_empty()
        True
        """
        with cfg.aux_solver as saux:
            return not saux.is_sat(self._make_f(cfg))

    def to_lp(self, vars):
        """
        Dump the polyhedron in .lp format.
        """
        def _to_lp1(x, t):
            """
            Dump one type (equation or inequalites) to string.
            """
            s = ""
            for i,h in enumerate(x):
                s += f"  {t}{i}:"
                for j,c in enumerate(h[1:]):
                    if c:
                        s += f" {'' if c < 0 else '+'}{c} {vars[j]}"
                s += f" {'>=' if t == 'ie' else '='} {-h[0]}\n"
            return s
        s = ("MAXIMIZE\nSubject To\n" +
            _to_lp1(self.eqns, "eq") +
            _to_lp1(self.ieqs, "ie") +
            "END\n")
        return s


multi_lp_sep = "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"

def list_to_lp(l, vars):
    """
    Dump a list of polyhedra as a list of .lp file contents, separated by multi_lp_sep.
    """
    if not l:
        s = "\\ no solution"
    else:
        s = ""
        for cnt,p in enumerate(l):
            s += f"\n{multi_lp_sep} file {cnt}\n\n"
            s += p.to_lp(vars)
        s += f"\n{multi_lp_sep} end\n\n"
    return s

def list_to_subsystems(l,vars):
    """
    Dump a list of polyhedra as a list of truncated subsystems, separated by multi_lp_sep.
    """
    if not l:
        s = "\\ no solution"
    else:
        s = ""
        for cnt,p in enumerate(l):
            s += f"\n{multi_lp_sep} file {cnt}\n\n"
            q=p.terms
            q.sort(key=lambda varn: varn[0][1])
            for i,elem in enumerate(q):
                if elem[0][1]==elem[1][1]:
                    s += f"{elem[0][0]} {elem[1][0]}\n"
                else:
                    s += f"|{elem[0][0]}| >= |{elem[1][0]}|\n"
        s += f"\n{multi_lp_sep} end\n\n"
    return s

# ------ load .lp files

def dim_expand(l, maxlen):
    """
    >>> dim_expand([[1],[2]], 2)
    [(1, 0), (2, 0)]
    >>> dim_expand([[]], 2)
    [(0, 0)]
    >>> dim_expand([[1],[2,2]], 2)
    [(1, 0), (2, 2)]
    """
    r = []
    for i in l:
        r.append(tuple(i + [0] * (maxlen - len(i))))
    return r


def read_lp_file(f):
    r"""
    Convert a .lp file to a list of input eqns and ieqs for phwrap.

    >>> read_lp_file([ "Subject To", "  ie0: -1 x1 >= -6", "  eq0: +1 x1 = -22", "BOUNDS" ])
    ([(22, 1)], [(6, -1)])
    >>> read_lp_file([ "Subject To", "  ie0: -1 x2 >= -6", "  ie0: +1 x1 +1 x3 >= 4", "  eq0: +1 x1 = -22", "  eq0: +1 x1 -1 x2 -2 x3 = +14", "BOUNDS" ])
    ([(22, 1, 0, 0), (-14, 1, -1, -2)], [(6, 0, -1, 0), (-4, 1, 0, 1)])
    """
    valid = False
    ign = True
    ie = []
    eq = []
    maxlen = 0
    for l in f:
        l = l.split("\\", 1)[0].strip()
        if l == "Subject To":
            ign = False
            valid = True
            continue
        if l == "BOUNDS":
            ign = True
            continue
        if l == "END":
            break
        if not l or ign:
            continue
        a = l.split(":")
        if len(a) != 2:
            print("Invalid constraint line: {}".format(l))
            continue
        a = a[1].split("=")
        if len(a) != 2:
            print("Invalid sense for constraint: {}".format(l))
            continue
        iseq = True
        if a[0][-1] == ">":
            iseq = False
            a[0] = a[0][:len(a[0])-1]                   # snip off last char
        vec = [-int(a[1])]      # r.h.s.
        a = a[0].split()
        for c,v in zip(a[::2], a[1::2]):
            if v[0] != "x":
                print("Invalid variable name: {}".format(v))
                break
            nb = int(v[1:])
            if len(vec) <= nb:
                vec = vec + [0] * (nb - len(vec) + 1)   # fill with zeros
            vec[nb] = int(c)
        if iseq:
            eq.append(vec)
        else:
            ie.append(vec)
        maxlen = max(maxlen, len(vec))

    if not valid:
        return None
    eq = dim_expand(eq, maxlen)
    ie = dim_expand(ie, maxlen)
    return eq, ie


def load_ph_from_lp(seq, l):
    """
    Read a sequence from an .lp file, create the polyhedron and append it to list 'l'.
    """
    r = read_lp_file(seq)
    if r is None:
        # not an .lp file
        return
    eqns, ieqs = r
    p = PhKeep(eqns=eqns, ieqs=ieqs)
    if p is not None:
        l.append(p)


def load_ph_from_multi_lp(fname):
    """
    Open and read (multi-) .lp file.
    Returns a list of polyhedra.
    """
    l = []
    with open(fname) as f:
        lines = []
        for ln in f:
            if ln.startswith(multi_lp_sep) and lines:
                # end of section, create polyhedron
                load_ph_from_lp(lines, l)
                lines = []
            else:
                # collect line
                lines.append(ln.rstrip())
        if lines:
            # seems like a non-multi .lp file, create polyhedron
            load_ph_from_lp(lines, l)
    return l


def tester():
    from .prt import prt
    prt("test: phkeep")
    import doctest, sys
    this_mod = sys.modules[__name__]    
    doctest.testmod(this_mod, verbose=False)
