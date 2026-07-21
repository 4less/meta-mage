#!/usr/bin/env python3
"""
Competitive rescue selection: re-admit guard-dropped markers that are CONTESTED.

The nucleotide guard (specificity_guard.py) drops a marker whenever its sequence
also occurs in an out-of-clade genome -- an ABSOLUTE presence test. nakedness.py
then asks the competitive question the guard never does: does every off-target
clade that carries the region ALSO carry its own marker there? If so the marker is
CONTESTED -- under a best-hit classifier the off-target reads map to the competing
marker, not this one, so the drop was conservative and the marker is safe to keep.

This selector pulls the CONTESTED markers' rows and sequences back out of the
PRE-guard emitted set (markers.emitted.tsv / markers.nuc.fasta -- the guard removed
them from markers.specific.*), so they can be merged onto the guard survivors.

NAKED markers (>= 1 off-target clade with no competitor -> genuine read theft) are
left dropped. Optionally restrict the rescue to one rank (default: all ranks).

Inputs
  --nakedness  nakedness.tsv  (cols: rank clade cluster ... verdict ...)
  --markers    markers.emitted.tsv  (FULL scored+emitted marker table, pre-guard)
  --fasta      markers.nuc.fasta    (FULL marker sequences, pre-guard)
Outputs
  --out_markers / --out_fasta  the contested subset (marker table + FASTA)
"""
import argparse


def load_contested_keys(path, rank=None):
    """(rank,clade,cluster) for every marker with verdict == 'contested'."""
    keys = set()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        vi = col.get("verdict")
        if vi is None:
            raise SystemExit("nakedness table lacks a 'verdict' column")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= vi or f[vi] != "contested":
                continue
            if rank and f[0] != rank:
                continue
            keys.add((f[0], f[1], f[2]))
    return keys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nakedness", required=True, help="nakedness.tsv (verdict col)")
    ap.add_argument("--markers", required=True,
                    help="pre-guard marker table (markers.emitted.tsv)")
    ap.add_argument("--fasta", required=True,
                    help="pre-guard marker FASTA (markers.nuc.fasta)")
    ap.add_argument("--rank", default="",
                    help="restrict rescue to this rank (default: all ranks)")
    ap.add_argument("--out_markers", required=True)
    ap.add_argument("--out_fasta", required=True)
    args = ap.parse_args()

    keys = load_contested_keys(args.nakedness, args.rank or None)

    # Pull the contested rows out of the pre-guard marker table. These markers were
    # dropped by the guard, so they are NOT in markers.specific.tsv -- markers.emitted
    # (EMIT_REPS output, before the guard) is where their rows still live.
    n_rows = 0
    with open(args.markers) as fh, open(args.out_markers, "w") as out:
        out.write(fh.readline())  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and (f[0], f[1], f[2]) in keys:
                out.write(line if line.endswith("\n") else line + "\n")
                n_rows += 1

    # Subset the pre-guard FASTA to the contested markers (headers are
    # rank|clade|cluster|genome_id; clades contain spaces, so key off the '|' split).
    written = 0
    with open(args.fasta) as fh, open(args.out_fasta, "w") as out:
        emit = False
        for line in fh:
            if line.startswith(">"):
                p = line[1:].rstrip("\n").split("|")
                key = (p[0], p[1], p[2]) if len(p) >= 3 else None
                emit = key in keys
                if emit:
                    out.write(line)
                    written += 1
            elif emit:
                out.write(line)

    print(f"select_contested: re-admitted {n_rows} contested markers "
          f"({written} sequences)"
          + (f" at rank {args.rank}" if args.rank else ""))


if __name__ == "__main__":
    main()
