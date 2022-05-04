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

import os, sys

from .util import mytime, used_mem, nice_combos, get_isodatetime, combos
from .phkeep import Bag, PhKeep, list_to_lp, list_to_subsystems
from .smt import MySolver, model_to_vec, Not, Or, And
from .preproc import preprocess
from .tropicalize import Tropicalize
from .config import Return


class SmtData:
    def __init__(self, *, space_dim):
        self.smt_time = 0
        self.isect_time = 0
        self.insert_time = 0
        self.space_dim = space_dim
        self.vars = [f"x{i+1}" for i in range(self.space_dim)]
        self.round = 1


def compute_dnf(cfg, bb, ret=None):
    """
    The main worker loop.
    cfg.minimize=True will minimize found polyhedra before inserting them into the result list.
    cfg.dump_ph=True will dump found points and polyhedra to the log file.
    """

    def status_print(len_lp, c=' '):
        # print status
        mem = used_mem()
        mem_str = "" if mem is None else f", mem {mem[0] >> 20}/{mem[1] >> 20}"
        t = mytime()
        cfg.prt(f"[{len_lp} ph{c}{mem_str}, smt {d.smt_time:.3f}, min {d.isect_time:.3f}, ins {d.insert_time:.3f}, tot {t-d.tt:.3f}]", flushfile=True, screen=cfg.verbose >= 1)

    d = SmtData(space_dim=bb[-1][0].space_dim())
    d.tt = mytime()
    d.status = status_print
    d.sub_solver = MySolver(name=cfg.solver_name, incremental=True)

    lp = Bag()
    lps = None

    with MySolver(name=cfg.solver_name, incremental=True) as solv:
        while True:
            cfg.prt(end=f"\nRunning SMT solver ({d.round})... ", screen=cfg.verbose >= 2)
            t = mytime()

            if lps:
                # exclude polyhedra from last pass
                ff = Not(cfg.mangler.smt_phl(lps))
            else:
                # first, add the whole system
                ff = cfg.mangler.smt_sys(bb)
            solv.add_assertion(ff)

            # solve!
            v = model_to_vec(solv.get_model(), d.space_dim) if solv.solve() else None
            t = mytime() - t
            d.smt_time += t
            cfg.prt(f"{t:.3f} sec", screen=cfg.verbose >= 2)

            if v is None:
                cfg.prt(f"No {'more ' if lp else ''}solutions found", screen=cfg.verbose >= 2 or not lp)
                break

            if cfg.dump_ph:
                cfg.prt(f"Found point ({','.join(str(i) for i in v)})", screen=False)

            # cycle through all bags to see which polyhedra contain the point
            cfg.prt(end=f"Checking... ", screen=cfg.verbose >= 2)
            t = mytime()
            one = PhKeep()
            bbs = []
            for b in bb:
                if len(b) == 1:
                    one.add_constraints(phkeep=b[0])
                    # if its just one polyhedron per bag, add constraints without checking for containment
                else:
                    # search which polyhedra contain the point
                    bs = Bag()
                    for p in b:
                        if p.contains_point(v):
                            bs.append(p)
                    if len(bs) == 1:
                        # add to general constraints if only one polyhedron contains the point
                        one.add_constraints(phkeep=bs[0])
                    else:
                        bbs.append(bs)

            # possibly minimize
            if cfg.minimize:
                cfg.prt(end=f"Minimizing ({len(one.eqns)},{len(one.ieqs)})... ", screen=cfg.verbose >= 2)
                one.minimize(cfg)
                cfg.prt(end=f"({len(one.eqns)},{len(one.ieqs)}), ", screen=cfg.verbose >= 2)
            t = mytime() - t
            d.isect_time += t
            cfg.prt(f"{t:.3f} sec", screen=cfg.verbose >= 2)

            if bbs:
                lps = compute_dnf_sub(cfg, d, bbs, v, len(lp), one)
            else:
                lps = [one]
                d.round += 1
                if cfg.dump_ph:
                    cfg.prt("\n\n" + one.to_lp(d.vars), screen=False)
                    #cfg.prt(f"{one._make_f(cfg).serialize()}", screen=False)
                d.status(len(lp) + 1)

            lp.extend(lps)

    if cfg.dump_ph:
        cfg.prt("\n\n\nSolution:\n" + list_to_lp(lp, d.vars), screen=False)
    cfg.prt(screen=cfg.verbose >= 1)

    #import collections
    #fvec = collections.Counter([p.dim() for p in lp])
    #fvec_str = str([fvec.get(i, 0) for i in range(max(fvec.keys())+1)])
    #if fvec:
    #    cfg.prt(f"f-vector: {fvec_str}", screen=cfg.verbose >= 1)
    total_time = mytime() - d.tt
    cfg.prt(f"{len(lp)} polyhedra, {d.round} rounds, smt {d.smt_time:.3f} sec, min {d.isect_time:.3f} sec, insert {d.insert_time:.3f}, total {total_time:.3f} sec", 
        flushfile=True, screen=cfg.verbose >= 1)

    if ret:
        ret.rounds = d.round
        #ret.fvector = fvec
        #ret.fvector_str = fvec_str
        ret.total_time = total_time
        ret.smt_time = d.smt_time
        ret.isect_time = d.isect_time
        ret.insert_time = d.insert_time
        ret.phs = len(lp)

    return lp


def compute_dnf_sub(cfg, d, bb, v, lp_len, one):
    round = 0
    lp, lp_max = Bag(), Bag()
    solv = d.sub_solver.solver
    solv.push()

    while True:
        if round == 0:
            t = mytime()
            solv.add_assertion(cfg.mangler.smt_sys(bb + [[one]]))
            t = mytime() - t
            d.smt_time += t
            t = mytime()
            two = PhKeep(phkeep=one)
            single = False
            # select some polyhedron from each bag
            for b in bb:
                two.add_constraints(phkeep=b[0])
            assert two.contains_point(v)
        else:
            cfg.prt(end=f"\nRunning SMT sub solver ({round})... ", screen=cfg.verbose >= 2)
            t = mytime()
            
            # exclude polyhedron from last pass
            solv.add_assertion(Not(cfg.mangler.smt_ph(two)))
            
            # solve!
            v = model_to_vec(solv.get_model(), d.space_dim) if solv.solve() else None
            t = mytime() - t
            d.smt_time += t
            cfg.prt(f"{t:.3f} sec", screen=cfg.verbose >= 2)

            if v is None:
                cfg.prt(f"No more solutions found", screen=cfg.verbose >= 2)
                break

            if cfg.dump_ph:
                cfg.prt(f"Found point ({','.join(str(i) for i in v)})", screen=False)

            # cycle through all bags, see which polyhedra contain the point
            cfg.prt(end=f"Checking... ", screen=cfg.verbose >= 2)
            t = mytime()
            two = PhKeep(phkeep=one)
            single = True
            for b in bb:
                # we know that all bags have at least 2 polyhedra
                assert len(b) > 1

                # search which polyhedra contain the point
                first = True
                for p in b:
                    if p.contains_point(v):
                        if first:
                            two.add_constraints(phkeep=p)
                            if not single:
                                break
                            first = False
                        else:
                            single = False
                            break

        if single and cfg.minimize:
            cfg.prt(end=f"Minimizing ({len(two.eqns)},{len(two.ieqs)})... ", screen=cfg.verbose >= 2)
            two.minimize(cfg, min_eqns=one.eqns)
            cfg.prt(end=f"({len(two.eqns)},{len(two.ieqs)}), ", screen=cfg.verbose >= 2)
        t = mytime() - t
        d.isect_time += t
        cfg.prt(f"{t:.3f} sec", screen=cfg.verbose >= 2)

        cfg.prt(end="Inserting... ", screen=cfg.verbose >= 2)
        if cfg.dump_ph:
            cfg.prt("\n\n" + two.to_lp(d.vars), screen=False)
            # cfg.prt(f"{two._make_f(cfg).serialize()}", screen=False)
        two.point = v
        t = mytime()
        if single:
            # if only one polyhedron contains this point, then no later polyhedra can contain this polyhedron.
            # otherwise, they would have contained this point as well.
            # still, earlier polyhedra could still be contained in this one.  hence, filter already found ones.
            lp.remove_included(cfg, two)
            lp_max.append(two)
        else:
            lp.insert_larger(cfg, two)
        t = mytime() - t
        d.insert_time += t
        cfg.prt(f"{t:.3f} sec", screen=cfg.verbose >= 2)

        d.status(lp_len + len(lp_max) + len(lp), '+' if round == 0 else '*')
        round += 1

    solv.pop()
    d.round += round
    
    for p in lp:
        if cfg.minimize:
            cfg.prt(end=f"Minimizing ({len(p.eqns)},{len(p.ieqs)})... ", screen=cfg.verbose >= 2)
            t = mytime()
            lp_max.append(p.minimize(cfg, min_eqns=one.eqns))
            t = mytime() - t
            d.isect_time += t
            cfg.prt(end=f"({len(p.eqns)},{len(p.ieqs)}), ", screen=cfg.verbose >= 2)
            cfg.prt(f"{t:.3f} sec", screen=cfg.verbose >= 2)
        else:
            lp_max.append(p)

    return lp_max


def smt_compare(cfg, bb1, bb2, ret=None):
    """
    Compare two conjunctions of disjunctions of polyhedra for set equality.
    """
    t = mytime()
    f1 = cfg.mangler.smt_sys(bb1)
    f2 = cfg.mangler.smt_sys(bb2)
    f = Or(And(f1, Not(f2)), And(f2, Not(f1)))
    with cfg.aux_solver as saux:
        res = not saux.is_sat(f)
    t = mytime() - t

    if ret:
        ret.compare_time = t
        ret.compare = res

    return res


def status_blob():
    """
    Output version numbers of Python and used libraries.
    """
    import platform
    s = ""
    s += f"Python {sys.version} {platform.machine()}{'' if __debug__ else ', (no debug)'}\n"
    try:
        import pysmt
    except ImportError:
        pass
    else:
        try:
            s += f"pysmt {pysmt.__version__}, "
        except AttributeError:
            s += f"pysmt {pysmt.git_version()}, "

    try:
        import gmpy2
    except ImportError:
        pass
    else:
        s += f"gmpy2 {gmpy2.version()}, "

    try:
        import sympy
    except ImportError:
        pass
    else:
        s += f"sympy {sympy.__version__}, "

    if s.endswith(", "):
        s = s[:-2]
    s += "\n"
    return s


def compute_dnf_from_input(cfg, basename, status):
    dir = os.path.join(cfg.dbdir, basename)         # directory where it happens

    # build flags string
    flags = ("Flags: " +
        f"eps={cfg.eps} " +
        f"p={cfg.p} " +
        f"solver={cfg.solver_name} " +
        "no-" * (not cfg.ppone) + "ppone " +
        "no-" * (not cfg.ppsmt) + "ppsmt " +
        "no-" * (not cfg.pp2x2) + "pp2x2 " +
        "no-" * (not cfg.dump_ph) + "dump " +
        "no-" * (not cfg.compare) + "comp " +
        f"maxadd={cfg.maxadd} " * cfg.pp2x2 +
        "")

    # build 'mainname' which reflecs all important settings
    mainname = f"ep{str(1/cfg.eps).replace('/','d')}"
    if cfg.p != 1:
        mainname += f"-p{cfg.p}"
    if cfg.bench_logs:
        nameadd = "o" * cfg.ppone + "s" * cfg.ppsmt + f"x{cfg.maxadd}" * cfg.pp2x2
        if nameadd:
            nameadd = "-" + nameadd
        import platform
        nameadd += f"-{cfg.solver_name}-{platform.system()}"
        mainname += nameadd

    if cfg.save_log:
        fname = os.path.join(dir, f"{mainname}-log.txt")
        cfg.prt.set_log_file(open(fname, "w"))
    else:
        cfg.prt.set_log_file()
    cfg.prt("\n----------------------------------------------------------------------", screen=cfg.verbose >= 1)
    cfg.prt(f"{status}\nDatetime: {get_isodatetime()}", screen=False)
    cfg.prt(f"Model: {basename}", screen=cfg.verbose >= 1)
    cfg.prt(f"Logfile: {fname if cfg.save_log else '-'}", screen=cfg.verbose >= 1)
    cfg.prt(f"{flags}\n", screen=cfg.verbose >= 1)

    # to return some data
    ret = Return()

    trop = Tropicalize(cfg)
    param_fn = os.path.join(dir, "parameters_q.txt")
    eqn_fns = [os.path.join(dir, "vector_field_q.txt"), os.path.join(dir, "conservation_laws_q.txt")]
    slow_fn = os.path.join(dir, "slow.txt")
    bb0 = bb = trop.tropicalize(param_fn, eqn_fns, slow_fn, ret)

    cfg.prt(f"Bag sizes: {' '.join(str(len(b)) for b in bb)}", screen=cfg.verbose >= 1)
    cfg.prt(f"Possible combinations: {nice_combos(bb)} in {len(bb)} bags", screen=cfg.verbose >= 1)
    ret.combos = combos(bb)

    # (maybe) do preprocessing
    t = mytime()
    bb = preprocess(cfg, bb)
    ret.preprocessing_time = mytime() - t
    cfg.prt.flush_all()

    if not bb:
        cfg.prt("Trivial: no solutions", screen=cfg.verbose >= 1)
        ret.solution = []
    elif len(bb) == 1 and len(bb[0]) == 1:
        cfg.prt("Trivial: one solution", screen=cfg.verbose >= 1)
        ret.solution = bb[0]
    else:
        # the big computation
        ret.solution = compute_dnf(cfg, bb, ret)

    ret.super_total_time = mytime() - t
    cfg.prt(f"Preprocessing {ret.preprocessing_time:.3f} sec, super total {ret.super_total_time:.3f} sec", flushfile=True, screen=cfg.verbose >= 1)
    
    # make sure they describe the same thing
    if cfg.compare:
        cfg.prt(end="Comparing... ", flush=True, screen=cfg.verbose >= 1)
        smt_compare(cfg, [ret.solution], bb0, ret)
        cfg.prt(end=f"{ret.compare_time:.3f} sec.  ", screen=cfg.verbose >= 1)
        if ret.compare:
            cfg.prt("Solution matches input.", screen=cfg.verbose >= 1)
        else:
            cfg.prt("\n*** SOLUTION DOES NOT MATCH INPUT! ***", screen=cfg.verbose >= 0)

    # save solution to file
    if cfg.save_solution:
        fname = os.path.join(cfg.dbdir, basename, f"{mainname}-solutions.txt")
        cfg.prt(f"\nSaving solution to {fname}", screen=cfg.verbose >= 1)
        s = list_to_lp(ret.solution, trop.vars)
        with open(fname, "w", newline="\n") as f:
            f.write(s)
        trackname = os.path.join(cfg.dbdir, basename, f"{mainname}-subsystems.txt")
        cfg.prt(f"\nSaving associated subsystems to {trackname}", screen=cfg.verbose >= 1)
        track = list_to_subsystems(ret.solution,trop.vars)
        with open(trackname, "w", newline="\n") as g:
            g.write(track)

    return ret
