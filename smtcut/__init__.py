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

from .config import Config
from .util import mytime, nice_combos, get_isodatetime
from .phkeep import PhKeep, list_to_lp, Bag, load_ph_from_lp, load_ph_from_multi_lp
from .prt import prt
from .tropicalize import Tropicalize
from .preproc import preprocess
from .main import compute_dnf, smt_compare, compute_dnf_from_input, status_blob


def tester():
    from .util import tester
    tester()
    from .smt import tester
    tester()
    from .phkeep import tester
    tester()
    from .tropicalize import tester
    tester()
