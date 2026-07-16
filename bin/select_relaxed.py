#!/usr/bin/env python3
"""
Adaptive rescue selection: add the fewest lower-prevalence markers needed.

For each flagged species, the extra candidates that PASSED the clean cross-map
guard (no off-target hit at all) are added in order of decreasing in-prevalence
(then score) until the species reaches the marker threshold or the candidates run
out. Higher-prevalence tiers are consumed first, so relaxation goes 0.7 -> 0.6 ->
0.5 only as far as needed. Cross-mapping markers are never here: the relaxed guard
kept only clean ones.

Emits the rescued marker rows (to append to the final marker table) and the
matching sequences subset from the clean relaxed FASTA.
"""
import argparse
from collections import defaultdict


def count_base(path, rank):
    n = defaultdict(int)
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0] == rank:
                n[f[1]] += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base_markers", required=True, help="final markers (0.8 run)")
    ap.add_argument("--extra_clean", required=True,
                    help="relaxed candidates that passed the clean guard")
    ap.add_argument("--extra_fasta", required=True, help="clean relaxed marker FASTA")
    ap.add_argument("--rank", default="species")
    ap.add_argument("--threshold", type=int, required=True)
    ap.add_argument("--out_markers", required=True)
    ap.add_argument("--out_fasta", required=True)
    args = ap.parse_args()

    base = count_base(args.base_markers, args.rank)

    # group extra clean candidates per species, carrying in_prev + score.
    by_sp = defaultdict(list)
    with open(args.extra_clean) as fh:
        header = fh.readline().rstrip("\n")
        cols = header.split("\t")
        ip_i = cols.index("in_prevalence")
        sc_i = cols.index("score") if "score" in cols else len(cols) - 1
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != args.rank:
                continue
            try:
                ip = float(f[ip_i])
                sc = float(f[sc_i])
            except (ValueError, IndexError):
                ip, sc = 0.0, 0.0
            by_sp[f[1]].append((ip, sc, line if line.endswith("\n") else line + "\n",
                                (f[0], f[1], f[2])))

    rescued_keys = set()
    added_per_sp = {}
    with open(args.out_markers, "w") as out:
        out.write(header + "\n")
        for sp, rows in by_sp.items():
            need = args.threshold - base.get(sp, 0)
            if need <= 0:
                continue
            rows.sort(key=lambda r: (-r[0], -r[1]))  # in_prev desc, then score
            added = 0
            for _ip, _sc, line, key in rows:
                if added >= need:
                    break
                out.write(line)
                rescued_keys.add(key)
                added += 1
            if added:
                added_per_sp[sp] = added

    # subset the clean relaxed FASTA to the rescued markers.
    written = 0
    with open(args.extra_fasta) as fh, open(args.out_fasta, "w") as out:
        emit = False
        for line in fh:
            if line.startswith(">"):
                p = line[1:].rstrip("\n").split("|")
                key = (p[0], p[1], p[2]) if len(p) >= 3 else None
                emit = key in rescued_keys
                if emit:
                    out.write(line)
                    written += 1
            elif emit:
                out.write(line)

    print(f"select_relaxed: rescued {len(rescued_keys)} markers across "
          f"{len(added_per_sp)} species ({written} sequences); "
          + ", ".join(f"{s}:+{n}" for s, n in sorted(added_per_sp.items())))


if __name__ == "__main__":
    main()
