#!/usr/bin/env python3
"""
Nakedness of cross-map-dropped markers -- competitive (best-hit) definition.

The nucleotide guard (specificity_guard.py) drops a marker whenever its sequence
also occurs in a genome OUTSIDE its clade -- an ABSOLUTE presence test. But under a
best-hit classifier a marker only actually steals reads when the off-target region
maps BETTER to it than to the off-target species' OWN representative marker. If the
off-target clade has a competing marker that matches the region at least as well,
its reads stay home: the cross-map is CONTESTED and the drop was conservative.

For each dropped marker A (clade Y) and each off-target clade X we compare two
identities over the shared region:

    id_off  = best identity of A into X's genomes         (crossmap.m8: A vs all_cds)
              -- how well X's reads would map to A
    id_comp = best identity of A to any marker of clade X  (self.m8: A vs markers)
              -- how well X's own marker represents the same region

X is CONTESTED (keeps its reads) when id_comp >= id_off; else X is NAKED (A wins
X's reads -> genuine theft). A marker is CONTESTED overall (safe to re-admit under a
competitive rule) when EVERY off-target clade is contested, else NAKED (correctly
dropped). Only pass==0 markers (those the guard dropped) are classified.

Note (proxy): id_comp uses marker-vs-marker identity id(A, M_X) as a stand-in for
"how well X represents the region." M_X is near-core in X, so this is a tight proxy;
it is the strongest comparison computable from the marker-vs-marker + marker-vs-
genome searches without simulating reads.

Inputs
  --self_hits / --selfmap    marker-vs-marker mmseqs m8 + m<N> id map
  --cross_hits / --crossmap  marker-vs-all_cds mmseqs m8 + m<N> id map
  --manifest                 idx -> lineage (target genome clade at the marker's rank)
  --specificity              specificity_report.tsv (pass column: which were dropped)
"""
import argparse
from collections import defaultdict

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_POS = {r: i for i, r in enumerate(RANKS)}


def load_idmap(path):
    """m<N> -> (rank, clade, cluster) from 'm<N> <TAB> rank|clade|cluster|genome'."""
    qmap = {}
    with open(path) as fh:
        for line in fh:
            qid, hdr = line.rstrip("\n").split("\t", 1)
            p = hdr.split("|")
            if len(p) >= 3:
                qmap[qid] = (p[0], p[1], p[2])
    return qmap


def load_manifest(path):
    """idx -> lineage tuple (one clade per rank)."""
    idx_lineage = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            idx_lineage[int(f[col["idx"]])] = tuple(f[col[r]] for r in RANKS)
    return idx_lineage


def norm_id(v):
    v = float(v)
    return v / 100.0 if v > 1.0 else v


def genome_idx(seq_id):
    return int(seq_id[1:].split("_", 1)[0])


def load_offtarget_ids(hits_path, qmap, idx_lineage, min_id, min_aln):
    """(rank,clade,cluster) -> {off_target_clade: best id_off}. id_off is how well
    the marker's nucleotide sequence matches genomes of the off-target clade."""
    off = defaultdict(lambda: defaultdict(float))
    with open(hits_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            q, t, pid, aln = f[0], f[1], f[2], f[3]
            ident = norm_id(pid)
            if ident < min_id or int(aln) < min_aln:
                continue
            key = qmap.get(q)
            if key is None:
                continue
            rank, clade, _ = key
            ri = RANK_POS.get(rank)
            if ri is None:
                continue
            try:
                lineage = idx_lineage.get(genome_idx(t))
            except (ValueError, IndexError):
                continue
            if lineage is None:
                continue
            tclade = lineage[ri]
            if tclade == clade or tclade == "NA":
                continue  # in-clade (incl. own copy)
            if ident > off[key][tclade]:
                off[key][tclade] = ident
    return off


def load_competitor_ids(hits_path, qmap, min_id, min_aln):
    """(rank,clade,cluster) -> {competing_clade: best id_comp}. id_comp is the best
    identity of a marker of that clade to this marker over the shared region."""
    comp = defaultdict(lambda: defaultdict(float))
    with open(hits_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            q, t, pid, aln = f[0], f[1], f[2], f[3]
            ident = norm_id(pid)
            if ident < min_id or int(aln) < min_aln:
                continue
            ka = qmap.get(q)
            kb = qmap.get(t)
            if ka is None or kb is None:
                continue
            if ka[0] != kb[0]:          # different rank
                continue
            if kb[1] == ka[1]:          # same clade (includes own copies)
                continue
            if ident > comp[ka][kb[1]]:
                comp[ka][kb[1]] = ident
    return comp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self_hits", required=True,
                    help="marker-vs-marker mmseqs m8: query target pident alnlen")
    ap.add_argument("--selfmap", required=True,
                    help="m<N> <TAB> rank|clade|cluster|genome for the self search")
    ap.add_argument("--cross_hits", required=True,
                    help="marker-vs-all_cds mmseqs m8 (crossmap.m8)")
    ap.add_argument("--crossmap", required=True,
                    help="m<N> <TAB> rank|clade|cluster|genome for the crossmap search")
    ap.add_argument("--manifest", required=True,
                    help="idx -> lineage (target genome clade lookup)")
    ap.add_argument("--specificity", required=True,
                    help="specificity_report.tsv (pass column marks guard-dropped)")
    ap.add_argument("--min_id", type=float, required=True)
    ap.add_argument("--min_aln", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    idx_lineage = load_manifest(args.manifest)
    cross_qmap = load_idmap(args.crossmap)
    self_qmap = load_idmap(args.selfmap)
    off_ids = load_offtarget_ids(args.cross_hits, cross_qmap, idx_lineage,
                                 args.min_id, args.min_aln)
    comp_ids = load_competitor_ids(args.self_hits, self_qmap,
                                   args.min_id, args.min_aln)

    with open(args.specificity) as fh, open(args.out, "w") as out:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        pass_i = col.get("pass", 3)
        out.write("rank\tclade\tcluster\tn_offtarget_clades\tn_contested\t"
                  "n_naked\tverdict\tnaked_clades\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= pass_i:
                continue
            rank, clade, cluster = f[0], f[1], f[2]
            if f[pass_i] != "0":        # only guard-dropped markers
                continue
            key = (rank, clade, cluster)
            offc = off_ids.get(key, {})
            if not offc:
                continue
            competitors = comp_ids.get(key, {})
            # A clade keeps its reads (contested) iff its competing marker matches
            # the region at least as well as this marker matches that clade.
            naked = sorted(x for x, id_off in offc.items()
                           if competitors.get(x, 0.0) < id_off)
            n_cont = len(offc) - len(naked)
            verdict = "naked" if naked else "contested"
            out.write(f"{rank}\t{clade}\t{cluster}\t{len(offc)}\t{n_cont}\t"
                      f"{len(naked)}\t{verdict}\t{';'.join(naked)}\n")


if __name__ == "__main__":
    main()
