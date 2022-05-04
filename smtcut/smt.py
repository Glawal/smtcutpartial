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

# is imported by: phkeep, config, main, preproc, smtcut

try:
    import pysmt
except ModuleNotFoundError:
    print(
"""
########################################
### Install pysmt: pip install pysmt ###
########################################
"""
    )
    raise

from pysmt.shortcuts import (Symbol, And, Or, Not, GE, Plus, Times, Equals, Real, Solver, FALSE)
from pysmt.typing import REAL
from pysmt.logics import QF_LRA
from pysmt.exceptions import NoSolverAvailableError

try:
    Solver(name="msat")
except NoSolverAvailableError:
    print(
"""
##########################################
### Install solvers with pysmt-install ###
##########################################
"""
    )
    raise


class MySolver:
    """
    Keep a solver instance that is not destroyed.  Just push() and pop() to keep assertions tidy.
    Works well with MSAT, Z3 and Yices.
    It is faster to use push()/pop() than to use reset_assertions().
    """
    def __init__(self, name="msat", **kwargs):
        self.solver = Solver(name, logic=QF_LRA, **kwargs)

    def __enter__(self, **kwargs):
        assert not kwargs
        self.solver.push()
        return self.solver

    def __exit__(self, type, value, traceback):
        self.solver.pop()


class FormulaMangler:
    def __init__(self):
        self.clear()

    def clear(self):
        self.eq_cache = {}
        self.ie_cache = {}

    def smt_point(self, v):
        """
        >>> m = FormulaMangler()
        >>> m.smt_point([1,2,-3])
        ((x1 = 1.0) & (x2 = 2.0) & (x3 = -3.0))
        """
        ll = []
        for cnt,val in enumerate(v, 1):
            ll.append(Equals(Symbol("x" + str(cnt), REAL), Real(val)))
        return And(ll)

    def smt_sum(self, d):
        """
        >>> m = FormulaMangler()
        >>> m.smt_sum([1,1,0])
        (x1, -1.0)
        >>> m.smt_sum([-1,1,1])
        ((x1 + x2), 1.0)
        >>> m.smt_sum([0,1,-1])
        (x1, x2)
        >>> m.smt_sum([0,2,0])
        ((x1 * 2.0), 0.0)
        >>> m.smt_sum([1,-2,-3])
        (1.0, ((x1 * 2.0) + (x2 * 3.0)))
        """
        lp, lm = [], []
        for cnt,val in enumerate(d[1:], 1):
            if val != 0:
                if val > 0:
                    ll = lp
                else:
                    ll = lm
                    val = -val
                if val == 1:
                    ll.append(Symbol("x" + str(cnt), REAL))
                else:
                    ll.append(Times(Symbol("x" + str(cnt), REAL), Real(val)))
        if d[0]:
            # for some strange reason, BM103 -ee -msat does much more rounds (3006 vs 2571)
            # if we use 'if d[0] < 0 here.  Seems like it doesn't like zero on the r.h.s. (BM103 -e11 -msat)
            if not lm:
                lm.append(Real(-d[0]))
            else:
                lp.append(Real(d[0]))
        return Plus(lp) if lp else Real(0), Plus(lm) if lm else Real(0)

    def smt_eqn(self, q):
        # caching of formulas yields some 6% speed-up (BM103 -e11 -msat)
        """
        >>> m = FormulaMangler()
        >>> m.smt_eqn((1,2,3))
        (((x1 * 2.0) + (x2 * 3.0)) = -1.0)
        >>> m.smt_eqn((1,2,-3))
        (((x1 * 2.0) + 1.0) = (x2 * 3.0))
        >>> m.smt_eqn((-1,2,3))
        (((x1 * 2.0) + (x2 * 3.0)) = 1.0)
        >>> m.smt_eqn((1,-2,-3))
        (1.0 = ((x1 * 2.0) + (x2 * 3.0)))
        """
        # sp, sm = smt_sum(q)
        # return Equals(sp, sm)

        try:
            return self.eq_cache[q]
        except KeyError:
            sp, sm = self.smt_sum(q)
            self.eq_cache[q] = c = Equals(sp, sm)
            return c

    def smt_ieq(self, q):
        """
        >>> m = FormulaMangler()
        >>> m.smt_ieq((1,2,3))
        (-1.0 <= ((x1 * 2.0) + (x2 * 3.0)))
        >>> m.smt_ieq((1,2,-3))
        ((x2 * 3.0) <= ((x1 * 2.0) + 1.0))
        >>> m.smt_ieq((-1,2,3))
        (1.0 <= ((x1 * 2.0) + (x2 * 3.0)))
        >>> m.smt_ieq((1,-2,-3))
        (((x1 * 2.0) + (x2 * 3.0)) <= 1.0)
        """
        # sp, sm = smt_sum(q)
        # return GE(sp, sm)

        try:
            return self.ie_cache[q]
        except KeyError:
            sp, sm = self.smt_sum(q)
            self.ie_cache[q] = c = GE(sp, sm)
            return c

    def smt_eqns(self, q):
        return And([self.smt_eqn(i) for i in q])

    def smt_ieqs(self, q):
        return And([self.smt_ieq(i) for i in q])

    def smt_ph2(self, eq, ie):
        """
        >>> m = FormulaMangler()
        >>> m.smt_ph2([(-1,1,0)], [(1,1,0), (2,0,2)])
        ((x1 = 1.0) & (-1.0 <= x1) & (-2.0 <= (2.0 * x2)))
        """
        ll = []
        for i in eq:
            ll.append(self.smt_eqn(i))
        for i in ie:
            ll.append(self.smt_ieq(i))
        return And(ll)

    def smt_ph(self, p):
        return self.smt_ph2(p.eqns, p.ieqs)

    def smt_phl(self, l):
        if not l:
            # this is false
            return FALSE()
        return Or([self.smt_ph(p) for p in l])

    def smt_sys(self, bb):
        return And([self.smt_phl(b) for b in bb])


def model_to_vec(model, dim):
    """
    Convert a model from smt solver to a list of numbers.
    """
    v = [0] * dim   # if model doesn't specify a variable it doesn't make a difference
    for s,i in model:
        s = s.symbol_name()
        i = i.constant_value()
        v[int(s[1:]) - 1] = i
    return tuple(v)                                     # faster and uses less memory


def tester():
    from .prt import prt
    prt("test: smt")
    import doctest, sys
    this_mod = sys.modules[__name__]    
    doctest.testmod(this_mod, verbose=False)
