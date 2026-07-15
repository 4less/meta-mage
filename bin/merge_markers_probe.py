#!/usr/bin/env python3
"""
Probe: does recursively merging a species with its most-overlapping neighbours
increase its marker count? Answered from existing counts.tsv alone -- no re-run.

Merging two species is just relabelling genomes' clade, so for any merged clade C:
    in_count_C  = sum of the members' per-species in_count (counts.tsv)
    size_C      = sum of the members' genome counts
and a gene is a marker for C when
    in_prev_C  = (size_C - (score_sum_C - in_score_C)) / size_C  >= min_in
    out_prev_C = (marker_total - in_count_C) / (N - size_C)      <= max_out
    size_C >= min_clade_size,  then capped at max_per_clade.
Merging RESCUES genes shared with the absorbed species (out_prev falls) but can
DEMOTE genes core to only one half (in_prev falls) -- so the marker count is
recomputed at every step and the whole trajectory is printed, peak marked.

At each step the "next overlapping species" is the not-yet-merged species in the
same genus that shares the most of the current clade's CORE genes (the genes whose
out_prev it inflates). Greedy, one merge at a time.

This is the amino-acid-cluster prevalence selection (what SCORE does); it does not
re-run the nucleotide cross-map guard -- but merging only ever relaxes that guard
too, so this is a lower bound on the real marker gain.
"""
import argparse
import os
import sys
from collections import defaultdict


def load_manifest(path, genus_filter):
    """species -> (size, score_sum, genus); N_total. score_sum=size if no
    completeness column (legacy manifest)."""
    size = defaultdict(int)
    score = defaultdict(float)
    sp_genus = {}
    n_total = 0
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        has_comp = "completeness" in col
        for line in fh:
            f = line.rstrip("\n").split("\t")
            n_total += 1
            sp = f[col["species"]]
            if sp == "NA":
                continue
            size[sp] += 1
            score[sp] += float(f[col["completeness"]]) if has_comp else 1.0
            sp_genus[sp] = f[col["genus"]]
    return size, score, sp_genus, n_total


def load_counts(path, genus_species):
    """cluster -> {species: (in_count, in_score)} and cluster -> marker_total,
    restricted to species-rank rows whose clade is in genus_species."""
    per = defaultdict(dict)
    total = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[1] != "species":
                continue
            clade = f[2]
            if clade not in genus_species:
                continue
            cluster = f[0]
            in_count = int(f[3])
            marker_total = int(f[4])
            in_score = float(f[5]) if len(f) > 5 else float(in_count)
            per[cluster][clade] = (in_count, in_score)
            total[cluster] = marker_total
    return per, total


def marker_set(clade_species, per, total, size, score, n_total, p):
    """Return the set of clusters that are markers for the merged clade."""
    size_C = sum(size[s] for s in clade_species)
    score_C = sum(score[s] for s in clade_species)
    if size_C < p["min_clade_size"]:
        return set(), set()
    out_denom = n_total - size_C
    if out_denom <= 0:
        return set(), set()
    markers = set()
    core = set()
    members = clade_species
    # Only clusters carried by >=1 member can be markers.
    for cluster, spmap in per.items():
        inc = isc = 0.0
        hit = False
        for s in members:
            v = spmap.get(s)
            if v:
                inc += v[0]
                isc += v[1]
                hit = True
        if not hit:
            continue
        in_prev = (size_C - (score_C - isc)) / size_C
        if in_prev > 1.0:
            in_prev = 1.0
        if in_prev >= p["min_in"]:
            core.add(cluster)
        out_prev = (total[cluster] - inc) / out_denom
        if in_prev >= p["min_in"] and out_prev <= p["max_out"]:
            markers.add(cluster)
    return markers, core


def next_neighbour(clade_species, core, per, size, candidates):
    """Species (not in clade) sharing the most of the clade's CORE clusters."""
    best, best_ov = None, 0
    for b in candidates:
        if b in clade_species:
            continue
        ov = sum(1 for cl in core if b in per.get(cl, ()))
        if ov > best_ov:
            best, best_ov = b, ov
    return best, best_ov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="nextflow outdir")
    ap.add_argument("--species", required=True)
    ap.add_argument("--genus", help="restrict merges to this genus "
                    "(default: the target species' genus)")
    ap.add_argument("--min_in", type=float, default=0.80)
    ap.add_argument("--max_out", type=float, default=0.02)
    ap.add_argument("--min_clade_size", type=int, default=3)
    ap.add_argument("--max_per_clade", type=int, default=100)
    ap.add_argument("--max_steps", type=int, default=10)
    args = ap.parse_args()

    p = dict(min_in=args.min_in, max_out=args.max_out,
             min_clade_size=args.min_clade_size, max_per_clade=args.max_per_clade)
    manifest = os.path.join(args.base, "manifest", "manifest.tsv")
    counts = os.path.join(args.base, "counts", "counts.tsv")

    size, score, sp_genus, n_total = load_manifest(manifest, None)
    if args.species not in size:
        sys.exit(f"{args.species}: not in manifest")
    genus = args.genus or sp_genus[args.species]
    genus_species = {s for s, g in sp_genus.items() if g == genus}
    sys.stderr.write(f"[merge] genus {genus}: {len(genus_species)} species; "
                     f"loading counts...\n")
    per, total = load_counts(counts, genus_species)
    candidates = genus_species

    def cap(n):
        return min(n, p["max_per_clade"])

    clade = [args.species]
    markers, core = marker_set(clade, per, total, size, score, n_total, p)
    print(f"{'step':>4}  {'merge (+species)':<40} {'genomes':>7} "
          f"{'markers':>7} {'capped':>6} {'Δ':>5}  rescued/demoted")
    print(f"{0:>4}  {args.species:<40} {sum(size[s] for s in clade):>7} "
          f"{len(markers):>7} {cap(len(markers)):>6} {'--':>5}")
    trajectory = [(tuple(clade), len(markers), cap(len(markers)))]

    prev_markers = markers
    for step in range(1, args.max_steps + 1):
        nb, ov = next_neighbour(clade, core, per, size, candidates)
        if not nb or ov == 0:
            sys.stderr.write("[merge] no further overlapping species.\n")
            break
        clade.append(nb)
        markers, core = marker_set(clade, per, total, size, score, n_total, p)
        rescued = len(markers - prev_markers)
        demoted = len(prev_markers - markers)
        delta = len(markers) - trajectory[-1][1]
        print(f"{step:>4}  +{nb:<39} {sum(size[s] for s in clade):>7} "
              f"{len(markers):>7} {cap(len(markers)):>6} {delta:>+5}  "
              f"+{rescued}/-{demoted}  (shared {ov} core genes)")
        trajectory.append((tuple(clade), len(markers), cap(len(markers))))
        prev_markers = markers

    peak = max(trajectory, key=lambda t: t[2])
    print(f"\nPeak (capped) markers = {peak[2]} at a {len(peak[0])}-species clade:")
    print("  " + "  +  ".join(peak[0]))
    base_c = trajectory[0][2]
    print(f"Baseline {args.species}: {base_c} capped markers "
          f"-> merged peak {peak[2]} ({peak[2]-base_c:+d}).")


if __name__ == "__main__":
    main()
