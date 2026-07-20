#!/usr/bin/env python3
"""
Core-gene rescue candidate selector (uniqueness dropped).

For FLAGGED low-marker species only, select EVERY core gene cluster
(in_prevalence >= min_in) regardless of out_prevalence -- i.e. drop gate (b).
The nucleotide cross-map guard WITH masking then becomes the sole specificity
test: a shared core gene is kept if a long enough clean (non-cross-mapping)
window survives (specificity_guard.py --target all_cds). This finds markers that
live in the *locally divergent* window of a gene that is otherwise shared across
the complex -- exactly the genes gate (b) throws away, because overall uniqueness
(out_prev) and local maskability are orthogonal.

Why restrict to flagged species: cross-mapping every core gene is expensive, and
only species that ended below the marker threshold need it (the same trigger the
relax and merge probes use).

Already-emitted candidates (markers.emitted.tsv) are excluded -- the main guard
already ruled on them; here we test only the core genes it never saw.

Output schema matches markers.tsv / relax_candidates.tsv so EMIT_REPS consumes it:
    rank  clade  cluster  in_prevalence  out_prevalence  in_count  clade_size  score
"""
import argparse
from collections import defaultdict


def load_clade_sizes(path):
    sizes = {}
    n_total = None
    with open(path) as fh:
        next(fh)
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


def load_flagged(path, rank):
    flagged = set()
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[col["flagged"]] == "yes":
                flagged.add(f[col["species"]])
    return flagged


def load_emitted(path):
    """(clade, cluster) already emitted+guarded in the main run -- skip these."""
    done = set()
    if not path:
        return done
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        ci, cl = col.get("clade"), col.get("cluster")
        if ci is None or cl is None:
            return done
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(ci, cl):
                done.add((f[ci], f[cl]))
    return done


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--clade_sizes", required=True)
    ap.add_argument("--low_marker", required=True, help="low_marker_species.tsv")
    ap.add_argument("--emitted", help="markers.emitted.tsv (exclude already guarded)")
    ap.add_argument("--rank", default="species")
    ap.add_argument("--min_in", type=float, default=0.80)
    ap.add_argument("--min_clade_size", type=int, default=3)
    ap.add_argument("--score_out_exp", type=float, default=1.0)
    ap.add_argument("--max_per_clade", type=int, default=0,
                    help="0 = keep every core gene; >0 caps per species by score")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sizes, n_total = load_clade_sizes(args.clade_sizes)
    flagged = load_flagged(args.low_marker, args.rank)
    emitted = load_emitted(args.emitted)

    by_clade = defaultdict(list)
    with open(args.counts) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[1] != args.rank:
                continue
            cluster, rank, clade = f[0], f[1], f[2]
            if clade not in flagged:
                continue
            if (clade, cluster) in emitted:
                continue
            entry = sizes.get((rank, clade))
            if entry is None or entry[0] < args.min_clade_size:
                continue
            size, score_sum = entry
            out_denom = n_total - size
            if out_denom <= 0:
                continue
            in_count = int(f[3]); marker_total = int(f[4]); in_score = float(f[5])
            in_prev = (size - (score_sum - in_score)) / size
            if in_prev > 1.0:
                in_prev = 1.0
            if in_prev < args.min_in:                 # CORE only -- no uniqueness gate
                continue
            out_prev = (marker_total - in_count) / out_denom
            # rank the most-unique-and-core first (more likely to survive masking),
            # but every core gene is kept unless a cap is set.
            score = in_prev * (1.0 - out_prev) ** args.score_out_exp
            by_clade[(rank, clade)].append(
                (score, cluster, in_prev, out_prev, in_count, size))

    kept = 0
    with open(args.out, "w") as out:
        out.write("rank\tclade\tcluster\tin_prevalence\tout_prevalence\t"
                  "in_count\tclade_size\tscore\n")
        for (rank, clade), rows in by_clade.items():
            rows.sort(reverse=True)
            sel = rows[: args.max_per_clade] if args.max_per_clade > 0 else rows
            for score, cluster, in_prev, out_prev, in_count, size in sel:
                out.write(f"{rank}\t{clade}\t{cluster}\t{in_prev:.4f}\t"
                          f"{out_prev:.4f}\t{in_count}\t{size}\t{score:.4f}\n")
                kept += 1
    print(f"select_core: {len(flagged)} flagged species, {kept} core-gene "
          f"candidates (uniqueness dropped) -> {args.out}")


if __name__ == "__main__":
    main()
