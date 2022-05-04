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

import itertools

from .prt import prt
from .phkeep import PhKeep, Bag
from .util import mytime, nice_combos
from .smt import *


def collect_onebags(cfg, bb, minimize=False):
    """
    Collect all bags with only one polyhedron.
    List of bags has a combined 1-bag, if existant, in front.
    """
    newl = []                                           # new list of bags
    one = PhKeep()                                      # common constraints

    for b in bb:
        if len(b) == 0:
            # one empty bag kills the whole solution
            bb = []
            return
        if len(b) == 1:
            # collect constraints
            one.add_constraints(phkeep=b[0])
        else:
            # copy whole bag
            newl.append(b)

    if one.is_universe():
        # no 1-bag found, nothing changed
        return 

    if minimize:
        cfg.prt(end=f"Minimizing ({len(one.eqns)},{len(one.ieqs)})... ", flush=True, screen=cfg.verbose >= 1)
        one.minimize(cfg)
        cfg.prt(f"({len(one.eqns)},{len(one.ieqs)})", screen=cfg.verbose >= 1)

    return [[one]] + newl


def preprocess_polyhedra_one(cfg, bb):
    """
    Do 1-bag preprocessing.
    """
    cfg.prt(f"\nFiltering bags...", screen=cfg.verbose >= 1)
    loop = True
    tt = mytime()

    # cycle as long as we have bags with one polyhedron
    while loop:
        cfg.prt("Collecting...", screen=cfg.verbose >= 1)
        bb = collect_onebags(cfg, bb)
        if not bb or len(bb[0]) > 1:
            break

        # check all polyhedra against those constraints
        loop = False
        assert len(bb[0]) == 1
        one = bb[0][0]
        bbnew = []
        for i, b in enumerate(bb[1:], 1):
            cfg.prt(end=f"Bag {i}/{len(bb)}, {len(b)} polyhedra", screen=cfg.verbose >= 1)
            bnew = Bag(name=b.name)      # new bag
            for p in b:
                if not p.test_intersection_empty(cfg, one):
                    bnew.append(p)
            if len(bnew) <= 1:
                # new bag has only one polyhedron: loop again
                loop = True
            bbnew.append(bnew)
            drop = len(b) - len(bnew)
            cfg.prt(f" => {drop} dropped" if drop else "", screen=cfg.verbose >= 1)
        bb = [[one]] + bbnew

    tt = mytime() - tt
    cfg.prt(f"Possible combinations after preprocessing: {nice_combos(bb)} in {len(bb)} bags", screen=cfg.verbose >= 1)
    cfg.prt(f"Bag sizes: {' '.join(str(len(b)) for b in bb)}", screen=cfg.verbose >= 1)
    cfg.prt(f"Time: {tt:.3f} sec", screen=cfg.verbose >= 1)


def preprocess_polyhedra_smt(cfg, bb):
    """
    Use SMT to check for superfluous polyhedra.
    Very effective, but might take long!
    """
    cfg.prt(f"\nChecking for superfluous polyhedra...", screen=cfg.verbose >= 1)

    tt = mytime()
    with cfg.aux_solver as saux:
        onecnt = 0
        # sort with largest bags in front
        bb = sorted(bb, key=lambda b: -len(b))
        for i in range(0, len(bb)):
            b = bb[i]       # work on this bag
            if not getattr(b,"ppsmt", False) and len(b) > 1:
                cfg.prt(f"Bag {i+1}/{len(bb)}, {len(b)} polyhedra:", screen=cfg.verbose >= 1)
                # init attribute to True
                for p in b:
                    p.keep = True
                    p.f = p._make_f(cfg)
                f_v = cfg.mangler.smt_sys(bb)     # the (current) prevariety

                # check all polyhedra
                for k,p in enumerate(b, 1):
                    # check if 'p' is required for 'bb', i.e. would 'bb' be the same if 'p' wasn't part of the set 'b'.
                    cfg.prt(end=f"  Polyhedron {k}... ", flush=True, screen=cfg.verbose >= 1)
                    t = mytime()
                    p.keep = False
                    A = Or([q.f for q in b if q.keep])  # b \ p
                    B = p.f
                    C = f_v
                    # does it hold that ((A \cup B) \cap C) == (A \cap C) ?
                    # since ((A \cup B) \cap C) is surely a superset of (A \cap C),
                    # we only need to test if ((A \cup B) \cap C) is a subset of )A \cap C).
                    # so we test if (((A \cup B) \cap C) \ (A \cap C)) is satisfiable.
                    # after some set operations that reduces to (B \cap C \cap Not(A)).
                    #f = And(Or([A, B]), C, Not(And(A, C)))
                    f = And([B, C, Not(A)])
                    #f = And([C, Not(A)])       # much slower: 386 vs 210

                    saux.push()
                    p.keep = saux.is_sat(f)
                    saux.pop()
                    t = mytime() - t
                    cfg.prt(f"{'required' if p.keep else 'superfluous'}, {t:.3f} sec", screen=cfg.verbose >= 1)

                # filter out non-required polyhedra
                bb[i] = Bag([p for p in b if p.keep])
                bb[i].ppsmt = True
                cfg.prt(f"  => {len(bb[i])} polyhedra left", screen=cfg.verbose >= 1)
                if not bb[i]:
                    return []       # bag empty
                # remove attribute
                for p in b:
                    del p.keep
                    p.f = None
            onecnt += len(bb[i]) == 1

    if onecnt > 1:
        bb = collect_onebags(cfg, bb)

    tt = mytime() - tt
    cfg.prt(f"Possible combinations after preprocessing: {nice_combos(bb)} in {len(bb)} bags", screen=cfg.verbose >= 1)
    cfg.prt(f"Bag sizes: {' '.join(str(len(b)) for b in bb)}", screen=cfg.verbose >= 1)
    cfg.prt(f"Time: {tt:.3f} sec", screen=cfg.verbose >= 1)

    return bb


class NothingDone(Exception):
    pass

def preprocess_polyhedra_2x2(cfg, bb):
    first = True
    tt = mytime()
    newl = []
    bb = sorted(bb, key=lambda b: len(b))

    # cycle through all elements from below
    while bb:
        b = bb[0]
        if len(b) < 2:  # too small
            newl.append(b)
            del bb[0]
            continue
        # first element found, search second one from above
        j = len(bb) - 1
        while j > 0:
            b1 = bb[j]
            if len(b1) * len(b) - len(b1) - len(b) <= cfg.maxadd:  # cheap enough
                break
            j -= 1
        else:
            if first:
                raise NothingDone
            break

        if first and cfg.verbose >= 1:
            cfg.prt(f"\nCombining small bags...")
        first = False
        del bb[j], bb[0]  # remove from old list
        # build new bag
        cfg.prt(end=f"Intersecting {len(b)}x{len(b1)}... ", screen=cfg.verbose >= 1)
        bnew = []
        for p,q in itertools.product(b, b1):
            c = p.intersect(q)
            if not c.test_empty(cfg):
                c.minimize(cfg)
                bnew.append(c)
        cfg.prt(f"{len(bnew)} polyhedra", screen=cfg.verbose >= 1)
        # re-insert
        bb.append(bnew)
        bb = sorted(bb, key=lambda b: len(b))

    newl.extend(bb)
    bb = collect_onebags(cfg, newl)

    tt = mytime() - tt
    cfg.prt(f"Possible combinations after preprocessing: {nice_combos(bb)} in {len(bb)} bags", screen=cfg.verbose >= 1)
    cfg.prt(f"Bag sizes: {' '.join(str(len(b)) for b in bb)}", screen=cfg.verbose >= 1)
    cfg.prt(f"Time: {tt:.3f} sec", screen=cfg.verbose >= 1)

    return bb


def preprocess(cfg, bb):
    """
    Do all the requested preprocessing.
    """
    if len(bb) > 1 and cfg.ppone:
        preprocess_polyhedra_one(cfg, bb)
    if len(bb) > 1 and cfg.ppsmt:
        bb = preprocess_polyhedra_smt(cfg, bb)
    if len(bb) > 1 and cfg.pp2x2:
        while True:
            try:
                bb = preprocess_polyhedra_2x2(cfg, bb)
            except NothingDone:
                break
            if len(bb) > 1 and cfg.ppsmt:
                bb = preprocess_polyhedra_smt(cfg, bb)
    if len(bb) > 1 and cfg.pponemin:
        cfg.prt(screen=cfg.verbose >= 1)
        bb = collect_onebags(cfg, bb, minimize=cfg.minimize)

    return bb
