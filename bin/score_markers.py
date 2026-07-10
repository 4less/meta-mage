#!/usr/bin/env python3
"""
Select clade-specific markers from per-clade counts.

For each (cluster, rank, clade):
    in_prevalence  = in_count / clade_size
    out_prevalence = (marker_total - in_count) / (N - clade_size)
Keep markers with in_prevalence >= min_in and out_prevalence <= max_out, on
clades of at least min_clade_size. Rank candidates per clade by (in - out) and
cap at max_per_clade to keep the classifier DB balanced.

Streaming: the counts file is read once; only the passing candidates are held
(far smaller than the full counts), then grouped and capped.
"""
import argparse
from collections import defaultdict


def load_clade_sizes(path):
    sizes = {}
    n_total = None
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            rank, clade, size = line.rstrip("\n").split("\t")
            if rank == "__TOTAL__":
                n_total = int(size)
            else:
                sizes[(rank, clade)] = int(size)
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
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sizes, n_total = load_clade_sizes(args.clade_sizes)

    # candidates grouped by (rank, clade) -> list of (score, row)
    by_clade = defaultdict(list)
    with open(args.counts) as fh:
        next(fh)  # header
        for line in fh:
            cluster, rank, clade, in_count, marker_total = line.rstrip("\n").split("\t")
            in_count = int(in_count)
            marker_total = int(marker_total)
            size = sizes.get((rank, clade))
            if size is None or size < args.min_clade_size:
                continue
            out_denom = n_total - size
            if out_denom <= 0:
                continue  # clade == whole dataset: no outside, specificity undefined
            in_prev = in_count / size
            out_prev = (marker_total - in_count) / out_denom
            if in_prev >= args.min_in and out_prev <= args.max_out:
                score = in_prev - out_prev
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
