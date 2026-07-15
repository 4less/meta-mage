#!/usr/bin/env python3
"""
Select clade-specific markers from per-clade counts.

For each (cluster, rank, clade):
    in_prevalence  = (clade_size - (score_sum - in_score)) / clade_size
    out_prevalence = (marker_total - in_count) / (N - clade_size)
in_prevalence is completeness-weighted: present genomes count fully, each MISSING
genome subtracts only its completeness (score_sum = sum of the clade's genome
completeness; in_score = sum over the genomes carrying this cluster). With no
completeness metadata every score is 1.0, so this reduces to in_count / clade_size.
Keep markers with in_prevalence >= min_in and out_prevalence <= max_out, on
clades of at least min_clade_size. Rank candidates per clade by
    score = in_prevalence * (1 - out_prevalence) ** score_out_exp
and cap at max_per_clade to keep the classifier DB balanced. The exponent
weights how hard off-target prevalence is penalised (1 = linear).

Streaming: the counts file is read once; only the passing candidates are held
(far smaller than the full counts), then grouped and capped.
"""
import argparse
from collections import defaultdict


def load_clade_sizes(path):
    """(rank, clade) -> (size, score_sum); plus N_total. score_sum defaults to
    size for legacy clade_sizes files that predate the completeness column."""
    sizes = {}
    n_total = None
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rank, clade, size = f[0], f[1], int(f[2])
            score_sum = float(f[3]) if len(f) > 3 else float(size)
            if rank == "__TOTAL__":
                n_total = size
            else:
                sizes[(rank, clade)] = (size, score_sum)
    if n_total is None:
        raise SystemExit("clade_sizes missing __TOTAL__ row")
    return sizes, n_total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--clade_sizes", required=True)
    ap.add_argument("--min_in", type=float, required=True)
    ap.add_argument("--max_out", type=float, required=True)
    ap.add_argument("--min_clade_size", type=int, required=True)
    ap.add_argument("--max_per_clade", type=int, required=True)
    ap.add_argument("--score_out_exp", type=float, default=1.0,
                    help="exponent on (1 - out_prevalence) in the rank score")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sizes, n_total = load_clade_sizes(args.clade_sizes)

    # candidates grouped by (rank, clade) -> list of (score, row)
    by_clade = defaultdict(list)
    with open(args.counts) as fh:
        next(fh)  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cluster, rank, clade = f[0], f[1], f[2]
            in_count = int(f[3])
            marker_total = int(f[4])
            in_score = float(f[5]) if len(f) > 5 else float(in_count)
            entry = sizes.get((rank, clade))
            if entry is None or entry[0] < args.min_clade_size:
                continue
            size, score_sum = entry
            out_denom = n_total - size
            if out_denom <= 0:
                continue  # clade == whole dataset: no outside, specificity undefined
            # Completeness-weighted: present count fully, absences discounted by
            # each missing genome's completeness (score_sum - in_score).
            in_prev = (size - (score_sum - in_score)) / size
            if in_prev > 1.0:
                in_prev = 1.0
            out_prev = (marker_total - in_count) / out_denom
            if in_prev >= args.min_in and out_prev <= args.max_out:
                score = in_prev * (1.0 - out_prev) ** args.score_out_exp
                by_clade[(rank, clade)].append(
                    (score, cluster, in_prev, out_prev, in_count, size)
                )

    with open(args.out, "w") as out:
        out.write("rank\tclade\tcluster\tin_prevalence\tout_prevalence\t"
                  "in_count\tclade_size\tscore\n")
        for (rank, clade), rows in by_clade.items():
            rows.sort(reverse=True)  # best score first
            for score, cluster, in_prev, out_prev, in_count, size in rows[: args.max_per_clade]:
                out.write(f"{rank}\t{clade}\t{cluster}\t{in_prev:.4f}\t"
                          f"{out_prev:.4f}\t{in_count}\t{size}\t{score:.4f}\n")


if __name__ == "__main__":
    main()
