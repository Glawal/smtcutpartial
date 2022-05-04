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

from sympy import Rational

from .prt import prt
from .smt import FormulaMangler, MySolver

class Config:
    def __init__(self):
        # default program settings

        # SMT settings
        self.solver_name = "msat"

        # preprocessing
        self.ppsmt = False
        self.pp2x2 = False
        self.ppone = False
        self.pponemin = True
        self.maxadd = 20
        self.minimize = True

        # tropicalization
        self.eps = Rational(1, 11)
        self.p = 1

        # output files
        self.save_solution = True
        self.save_log = True
        self.bench_logs = False
        self.dump_ph = False

        # other options
        self.verbose = 1
        self.compare = False
        self.dbdir = ""

        # globally used objects
        self.prt = prt
        self.mangler = FormulaMangler()
        self.aux_solver = MySolver(name=self.solver_name, generate_models=False)


class Return:
    def __init__(self):
        self.compare = None
        self.tropicalization_time = None
        self.super_total_time = None
        self.phs = None
        self.preprocessing_time = None
        self.rounds = None
        self.combos = None
        self.total_time = None
        self.insert_time = None
        self.isect_time = None
        self.smt_time = None
        self.compare_time = None
        self.solution = None
