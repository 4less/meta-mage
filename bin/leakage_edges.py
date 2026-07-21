#!/usr/bin/env python3
"""
Turn a run's crossmap.m8 into a directed species->species leakage edge list for
leakage_tree.py. A "leak" is a marker gene of the query species whose CDS
cross-maps (>= min_id, >= min_aln) onto a genome belonging to a DIFFERENT species.

Resolves ids exactly as the pipeline does:
  * query m-id  -> query species   via query.map  (field 2: rank|clade|S:seq|acc)
  * target g<idx>_... -> species    via manifest idx column

Output (TSV): query_species  target_species  genes  hits
  genes = distinct query marker genes leaking into the target species (this is the
          weight leakage_tree.py draws); hits = total qualifying alignment rows.
"""
import argparse
from collections import defaultdict


def load_manifest(path):
    idx2sp = {}
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sp = f[col["species"]]
            if sp.startswith("s__"):
                sp = sp[3:]
            idx2sp[int(f[col["idx"]])] = sp
    return idx2sp


def load_query_species(path):
    m2sp = {}
    with open(path) as fh:
        for line in fh:
            mid, header = line.rstrip("\n").split("\t", 1)
            parts = header.split("|")
            if len(parts) >= 2:
                sp = parts[1]
                if sp.startswith("s__"):
                    sp = sp[3:]
                m2sp[mid] = sp
    return m2sp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crossmap", required=True)
    ap.add_argument("--query_map", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--min_id", type=float, default=0.95)
    ap.add_argument("--min_aln", type=int, default=100)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    idx2sp = load_manifest(args.manifest)
    m2sp = load_query_species(args.query_map)
    min_pid = args.min_id * 100.0

    genes = defaultdict(set)   # (qsp, tsp) -> {mid}
    hits = defaultdict(int)    # (qsp, tsp) -> count
    seen = skipped = 0
    with open(args.crossmap) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            mid, tgt = f[0], f[1]
            try:
                pid = float(f[2]); aln = int(f[3])
            except ValueError:
                continue
            if pid < min_pid or aln < args.min_aln:
                continue
            qsp = m2sp.get(mid)
            if qsp is None:
                skipped += 1
                continue
            try:
                tidx = int(tgt[1:].split("_", 1)[0])
            except (ValueError, IndexError):
                continue
            tsp = idx2sp.get(tidx)
            if tsp is None or tsp == qsp:      # unresolved or on-target (self)
                continue
            key = (qsp, tsp)
            genes[key].add(mid)
            hits[key] += 1
            seen += 1

    rows = sorted(genes, key=lambda k: (-len(genes[k]), -hits[k]))
    with open(args.out, "w") as out:
        out.write("query_species\ttarget_species\tgenes\thits\n")
        for k in rows:
            out.write(f"{k[0]}\t{k[1]}\t{len(genes[k])}\t{hits[k]}\n")

    print(f"leakage_edges: {seen} off-target hits over {len(rows)} directed "
          f"species pairs ({skipped} rows with unmapped query id) -> {args.out}")


if __name__ == "__main__":
    main()
