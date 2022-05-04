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

import sys
import os
from sympy import Rational

from smtcut import Config, prt, tester, compute_dnf_from_input, status_blob


def proc_options(argv, cfg, names=None):
    """
    Process the command line.
    """
    collect = 0

    for i in argv:
        if collect == 0:
            if i.startswith("-"):
                if i == "-h":
                    cfg.help = 1
                elif i == "--hh":
                    cfg.help = 2
                elif i == "--no-ppsmt":
                    cfg.ppsmt = False
                elif i == "--ppsmt":
                    cfg.ppsmt = True
                elif i == "--no-pp2x2":
                    cfg.pp2x2 = False
                elif i == "--pp2x2":
                    cfg.pp2x2 = True
                elif i == "--no-ppone":
                    cfg.ppone = False
                elif i == "--ppone":
                    cfg.ppone = True
                elif i == "--no-pponemin" or i == "--no-ppminone":
                    cfg.pponemin = False
                elif i == "--pponemin" or i == "--ppminone":
                    cfg.pponemin = True
                elif i == "--no-pp":
                    cfg.ppsmt = False
                    cfg.pp2x2 = False
                    cfg.ppone = False
                    cfg.pponemin = False
                elif i == "--pp":
                    cfg.ppsmt = True
                    cfg.pp2x2 = True
                    cfg.ppone = True
                    cfg.pponemin = True
                elif i == "--min":
                    cfg.minimize = True
                elif i == "--nomin" or i == "--no-min":
                    cfg.minimize = False
                elif i == "--nodump" or i == "--no-dump":
                    cfg.dump_ph = False
                elif i == "--dump":
                    cfg.dump_ph = True
                elif i == "--nocomp" or i == "--no-comp":
                    cfg.compare = False
                elif i == "--comp":
                    cfg.compare= True
                elif i == "--maxadd":
                    collect = 1
                elif i == "--smt":
                    collect = 2
                elif i in ("--msat", "--z3", "--cvc4", "--yices"):
                    cfg.solver_name = i[2:]
                elif i == "-d":
                    collect = 3
                elif i ==  "-e":
                    collect = 4
                elif i.startswith("-e"):
                    cfg.eps = Rational(1, i[2:])
                elif i ==  "-p":
                    collect = 5
                elif i.startswith("-p"):
                    cfg.p = int(i[2:])
                elif i.startswith("-q"):
                    cfg.verbose = 0
                elif i.startswith("-v"):
                    cfg.verbose += 1
                    for c in i[2:]:
                        if c == "0":
                            cfg.verbose = 0
                        elif c == "v":
                            cfg.verbose += 1
                elif i == "--bench-logs":
                    cfg.bench_logs = True
                elif i == "--no-bench-logs":
                    cfg.bench_logs = False
                elif i == "--save":
                    cfg.save_solution = True
                elif i == "--no-save":
                    cfg.save_solution = False
                elif i == "--log":
                    cfg.save_log = True
                elif i == "--no-log":
                    cfg.save_log = False
                elif i in ("-T", "--test"):
                    cfg.test = True
                elif i == "--no-test":
                    cfg.test = False
                else:
                    prt(f"*** Unrecognized option: {i}")
                    sys.exit(1)
            else:
                names.append(i)
        else:
            if collect == 1:
                cfg.maxadd = int(i)
            elif collect == 2:
                cfg.solver_name = i
            elif collect == 3:
                cfg.dbdir = i
            elif collect == 4:
                cfg.eps = Rational(1, i)
            elif collect == 5:
                cfg.p = int(i)
            collect = 0


def main():
    hello = f"This is SMTcut v5.3.3 by Christoph Lueders -- http://wrogn.com\n\n"
    status = hello + status_blob()

    cfg = Config()
    cfg.help = 0
    cfg.test = False

    names = []
    try:
        import shlex
        argv = shlex.split(os.environ["SMTCUT_OPTIONS"])
    except KeyError:
        argv = []
    argv.extend(sys.argv[1:])
    proc_options(argv, cfg, names)

    if cfg.help or (not names and not cfg.test):
        prt(end=hello +
f"""
Usage: {os.path.basename(sys.argv[0])} <model>... [<option>...]

General options (default: *):
  --(no-)dump     (*don't) dump results to log file
  --(no-)comp     (*don't) compare solution with input

  -e EPS          tropicalize with epsilon=1/EPS (default: 11)
  -p P            round with round(x*P)/P  (default: 1)
  -d DIR          set base directory to DIR (default: '')
  -q              be quiet
  -v              increase verbosity
  --smt NAME      use SMT solver 'NAME' (default: msat)
  --hh            more help
"""[1:] +
"""
Preprocessing options:
  --(no-)ppone    (*don't) preprocess 1-bags and intersect with other bags
  --(no-)ppsmt    (*don't) preprocess to remove unneeded polyhedra
  --(no-)pp2x2    (*don't) preprocess to combine small bags
  --(no-)pponemin (don't) *minimize 1-bag after preprocess
  --(no-)pp       *don't preprocess at all / do all preprocessing
  --(no-)min      (don't) *minimize found polyhedra

  --maxadd MAX    allow <= MAX new ph's for pp2x2 (default: 20)

Esoteric options:
  --(no-)save     (don't) *save solution
  --(no-)bench-logs  (*don't) use more specific output filenames for benchmarking
""" * (cfg.help > 1))
        sys.exit(1)

    prt(end=status, screen=cfg.verbose >= 1)

    if cfg.test:
        prt()
        tester()
        prt()

    # make sure spaces in options survive
    import re
    qc = ' '.join(f'"{i}"' if re.search(r'\s', i) else i for i in argv)
    status += f"Command line: {qc}\n"

    for basename in names:
        compute_dnf_from_input(cfg, basename, status)


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=False)

    main()
