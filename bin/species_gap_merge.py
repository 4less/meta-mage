#!/usr/bin/env python3
"""
Agglomerative "merge until the ANI gap turns positive" probe.

Starting from a seed species, repeatedly absorb the single nearest remaining
species (highest between-ANI) and recompute the gap

    gap(C) = within_low_quantile(C) - max over external E of between_high_quantile(C, E)

until either gap > 0 (C is a clean sequence-discrete cluster -- an SGB with a real
valley to everything else) or every species has been absorbed (no positive gap
exists: the taxon is a continuum, and there is no data-driven boundary, only a
choice of resolution).

Extreme min/max are outlier-dominated at genome scale, so the gap uses quantiles:
--within_q (default 0.05) for within-cohesion and --between_q (default 0.95) for
the nearest-neighbour similarity. If --counts is given, the marker count of the
growing clade is printed alongside (prevalence stage, from counts.tsv) so you can
watch resolution trade against marker payoff.

Input ANI is a skani -E sparse table (Ref Query ANI af_ref af_query); produce it
with the pipeline's ANI_SKANI step or `skani triangle -E` over the genus genomes.
"""
import argparse
import os
import statistics
import sys
from collections import defaultdict


def load_manifest(path):
    gid_sp, name2gid = {}, {}
    sp_genus = {}
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sp, genus = f[col["species"]], f[col["genus"]]
            if sp == "NA":
                continue
            gid = f[col["genome_id"]]
            gid_sp[gid] = sp
            sp_genus[sp] = genus
            name2gid[os.path.basename(f[col["path"]])] = gid
            name2gid[f[col["path"]]] = gid
    return gid_sp, sp_genus, name2gid


def parse_sparse(path, name2gid, gid_sp):
    """(species_a, species_b) sorted -> list of genome-pair ANIs (min_af filtered
    upstream by skani)."""
    pair = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("Ref_file", "#")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            a = name2gid.get(f[0]) or name2gid.get(os.path.basename(f[0]))
            b = name2gid.get(f[1]) or name2gid.get(os.path.basename(f[1]))
            if not a or not b:
                continue
            sa, sb = gid_sp.get(a), gid_sp.get(b)
            if not sa or not sb:
                continue
            try:
                ani = float(f[2])
            except ValueError:
                continue
            pair[tuple(sorted((sa, sb)))].append(ani)
    return pair


def quant(vals, q):
    if not vals:
        return None
    s = sorted(vals)
    if len(s) == 1:
        return s[0]
    i = q * (len(s) - 1)
    lo = int(i)
    frac = i - lo
    hi = min(lo + 1, len(s) - 1)
    return s[lo] * (1 - frac) + s[hi] * frac


def within_low(clade, pair, q):
    vals = []
    members = list(clade)
    for i, a in enumerate(members):
        for b in members[i:]:
            vals += pair.get(tuple(sorted((a, b))), [])
    return quant(vals, q)


def between_high(clade, ext, pair, q):
    vals = []
    for a in clade:
        vals += pair.get(tuple(sorted((a, ext))), [])
    return quant(vals, q)


# ------- optional marker count (prevalence stage), reused idea from the probe ---
def load_counts_species(path, species_set):
    per = defaultdict(dict)
    total = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[1] != "species" or f[2] not in species_set:
                continue
            per[f[0]][f[2]] = (int(f[3]), float(f[5]) if len(f) > 5 else float(f[3]))
            total[f[0]] = int(f[4])
    return per, total


def marker_count(clade, per, total, size, score, n_total, p):
    size_C = sum(size[s] for s in clade)
    score_C = sum(score[s] for s in clade)
    if size_C < p["min_clade_size"] or n_total - size_C <= 0:
        return 0
    n = 0
    for cl, spmap in per.items():
        inc = isc = 0.0
        hit = False
        for s in clade:
            v = spmap.get(s)
            if v:
                inc += v[0]; isc += v[1]; hit = True
        if not hit:
            continue
        in_prev = min(1.0, (size_C - (score_C - isc)) / size_C)
        out_prev = (total[cl] - inc) / (n_total - size_C)
        if in_prev >= p["min_in"] and out_prev <= p["max_out"]:
            n += 1
    return min(n, p["max_per_clade"])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--sparse", required=True, help="skani -E sparse ANI table")
    ap.add_argument("--species", required=True, help="seed species")
    ap.add_argument("--within_q", type=float, default=0.05)
    ap.add_argument("--between_q", type=float, default=0.95)
    ap.add_argument("--counts", help="counts.tsv -> also print marker count")
    ap.add_argument("--min_in", type=float, default=0.80)
    ap.add_argument("--max_out", type=float, default=0.02)
    ap.add_argument("--min_clade_size", type=int, default=3)
    ap.add_argument("--max_per_clade", type=int, default=100)
    ap.add_argument("--max_steps", type=int, default=50)
    args = ap.parse_args()

    gid_sp, sp_genus, name2gid = load_manifest(args.manifest)
    if args.species not in sp_genus:
        sys.exit(f"{args.species}: not in manifest")
    genus = sp_genus[args.species]
    genus_species = {s for s, g in sp_genus.items() if g == genus}
    pair = parse_sparse(args.sparse, name2gid, gid_sp)

    per = total = size = score = None
    p = dict(min_in=args.min_in, max_out=args.max_out,
             min_clade_size=args.min_clade_size, max_per_clade=args.max_per_clade)
    n_total = sum(1 for _ in gid_sp)
    if args.counts:
        size = defaultdict(int); score = defaultdict(float)
        for _g, s in gid_sp.items():
            size[s] += 1; score[s] += 1.0     # v1 manifest: unweighted
        per, total = load_counts_species(args.counts, genus_species)

    clade = {args.species}
    print(f"seed {args.species}  (genus {genus}, {len(genus_species)} species)")
    print(f"{'step':>4}  {'+species':<40} {'within':>7} {'nearest':>28} "
          f"{'betw':>6} {'gap':>7} {'markers':>7}")
    for step in range(args.max_steps + 1):
        wl = within_low(clade, pair, args.within_q)
        # nearest external species
        nearest, bh = None, -1e9
        for e in genus_species - clade:
            v = between_high(clade, e, pair, args.between_q)
            if v is not None and v > bh:
                nearest, bh = e, v
        gap = (wl - bh) if (wl is not None and nearest is not None) else float("nan")
        mk = (marker_count(clade, per, total, size, score, n_total, p)
              if args.counts else "-")
        added = "" if step == 0 else args._last
        near_s = (nearest or "-").replace("s__", "")[:26]
        wl_s = f"{wl:.2f}" if wl is not None else "NA"
        gap_s = f"{gap:+.2f}" if gap == gap else "NA"
        print(f"{step:>4}  {('+'+added if added else args.species):<40} "
              f"{wl_s:>7} {near_s:>28} {bh:>6.2f} {gap_s:>7} {str(mk):>7}")
        if wl is not None and nearest is not None and gap > 0:
            print(f"\n=> POSITIVE GAP reached: {args.species} resolves as a "
                  f"{len(clade)}-species cluster with a {gap:+.2f} ANI valley to "
                  f"the rest.")
            return
        if not nearest:
            break
        clade.add(nearest)
        args._last = nearest
    print(f"\n=> NO positive gap: merging ran to the whole genus "
          f"({len(clade)} species) without a valley opening. This taxon is an ANI "
          f"continuum -- there is no data-driven boundary, only a chosen resolution.")


if __name__ == "__main__":
    main()
