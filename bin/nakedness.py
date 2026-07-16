#!/usr/bin/env python3
"""
Nakedness of cross-map-dropped markers.

The nucleotide guard (specificity_guard.py) drops a marker whenever its sequence
also occurs in a genome OUTSIDE its clade. But a dropped marker is only genuinely
dangerous -- i.e. only actually steals reads and inflates the wrong taxon -- when
the off-target region is NOT itself covered by a marker of the off-target clade.
If the off-target clade DOES have a marker that matches (reads from that clade map
to both markers, so they compete/are shared), the cross-map is CONTESTED and the
drop was conservative: under a competitive classifier the marker could be kept.

For each off-target clade of a dropped marker we ask: does that clade have a
marker whose nucleotide sequence matches this marker (>= min_id over >= min_aln)?
That is answered by a marker-vs-marker self search (self_hits). Then:

    offtarget_clades  -- clades whose GENOMES carry the region (specificity_report)
    contested_clades  -- off-target clades that have a competing MARKER (self_hits)
    naked_clades      = offtarget_clades - contested_clades

A marker is CONTESTED (safe to keep under a competitive rule) when naked_clades is
empty, else NAKED (correctly dropped: >= 1 off-target clade steals reads with no
competitor). Only pass==0 markers (those the guard dropped) are classified.
"""
import argparse
import os
from collections import defaultdict


def load_selfmap(path):
    """m<N> -> (rank, clade, cluster) from 'm<N> <TAB> rank|clade|cluster|genome'."""
    qmap = {}
    with open(path) as fh:
        for line in fh:
            qid, hdr = line.rstrip("\n").split("\t", 1)
            p = hdr.split("|")
            if len(p) >= 3:
                qmap[qid] = (p[0], p[1], p[2])
    return qmap


def norm_id(v):
    v = float(v)
    return v / 100.0 if v > 1.0 else v


def load_contested(hits_path, qmap, min_id, min_aln):
    """(rank,clade,cluster) -> {off-target clade that has a competing marker}.

    A self hit query=marker A, target=marker B at the same rank but a different
    clade means clade(B) markers the region A cross-maps to, so A is contested by
    clade(B). In-clade hits (B in A's clade, incl. A's own copies) are ignored."""
    contested = defaultdict(set)
    with open(hits_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            q, t, pid, aln = f[0], f[1], f[2], f[3]
            if norm_id(pid) < min_id or int(aln) < min_aln:
                continue
            ka = qmap.get(q)
            kb = qmap.get(t)
            if ka is None or kb is None:
                continue
            if ka[0] != kb[0]:          # different rank
                continue
            if kb[1] == ka[1]:          # same clade (includes own copy)
                continue
            contested[ka].add(kb[1])    # clade(B) contests marker A
    return contested


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self_hits", required=True,
                    help="marker-vs-marker mmseqs m8: query target pident alnlen")
    ap.add_argument("--selfmap", required=True,
                    help="m<N> <TAB> rank|clade|cluster|genome for the marker fasta")
    ap.add_argument("--specificity", required=True,
                    help="specificity_report.tsv (needs the offtarget_clades column)")
    ap.add_argument("--min_id", type=float, required=True)
    ap.add_argument("--min_aln", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    qmap = load_selfmap(args.selfmap)
    contested = load_contested(args.self_hits, qmap, args.min_id, args.min_aln)

    with open(args.specificity) as fh, open(args.out, "w") as out:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        oc_i = col.get("offtarget_clades")
        if oc_i is None:
            raise SystemExit("specificity_report lacks offtarget_clades column "
                             "(rerun the guard); cannot assess nakedness")
        out.write("rank\tclade\tcluster\tn_offtarget_clades\tn_contested\t"
                  "n_naked\tverdict\tnaked_clades\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= oc_i:
                continue
            rank, clade, cluster, passed = f[0], f[1], f[2], f[3]
            if passed != "0":           # only guard-dropped markers
                continue
            off = set(x for x in f[oc_i].split(";") if x)
            if not off:
                continue
            cont = contested.get((rank, clade, cluster), set())
            naked = off - cont
            n_cont = len(off & cont)
            verdict = "contested" if not naked else "naked"
            out.write(f"{rank}\t{clade}\t{cluster}\t{len(off)}\t{n_cont}\t"
                      f"{len(naked)}\t{verdict}\t{';'.join(sorted(naked))}\n")


if __name__ == "__main__":
    main()
