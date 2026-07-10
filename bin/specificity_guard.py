#!/usr/bin/env python3
"""
Nucleotide specificity guard -- MetaPhlAn's final marker uniqueness check.

A marker earns its clade label from PROTEIN clustering, but reads are assigned in
NUCLEOTIDE space. A marker is only trustworthy if its nucleotide sequence is NOT
also present in a genome outside its target clade: otherwise reads from that
off-target taxon map to the marker and inflate the clade's abundance.

We align every emitted marker sequence against the species-rep CDS universe
(reps.ffn) and inspect the hits:
  * A hit to a genome whose clade AT THE MARKER'S RANK equals the marker's clade
    is in-clade (this includes the marker's own source copy) -> ignored.
  * A hit to a genome in a DIFFERENT clade, at >= identity over >= a read-length
    window, is an off-target cross-map -> the whole marker (rank,clade,cluster)
    is dropped.

Inputs are the mmseqs hit table (query, target, pident, alnlen), a query id map
(mmseqs truncates ids at whitespace, so queries are re-keyed to m<N> upstream),
the scored marker table, and the marker FASTA. Outputs are the marker table and
FASTA filtered to the surviving markers, plus a per-marker report.

Query id map:  m<N> <TAB> rank|clade|cluster|genome_id   (cluster keeps L:/S: tag)
Target header: g<idx>_<original CDS id>
"""
import argparse

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_POS = {r: i for i, r in enumerate(RANKS)}


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


def load_idmap(path):
    """m<N> -> (rank, clade, cluster) parsed from the original query header."""
    qmap = {}
    with open(path) as fh:
        for line in fh:
            qid, header = line.rstrip("\n").split("\t", 1)
            parts = header.split("|")
            if len(parts) >= 3:
                qmap[qid] = (parts[0], parts[1], parts[2])
    return qmap


def genome_idx(seq_id):
    # g<idx>_<original token>
    return int(seq_id[1:].split("_", 1)[0])


def norm_identity(value):
    """mmseqs may report pident as a fraction (0-1) or a percentage (0-100)."""
    v = float(value)
    return v / 100.0 if v > 1.0 else v


def scan_hits(path, qmap, idx_lineage, min_id, min_aln):
    """Return {(rank,clade,cluster): (n_offtarget, max_id, worst_clade)}."""
    offenders = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            query, target, pident, alnlen = f[0], f[1], f[2], f[3]
            if norm_identity(pident) < min_id or int(alnlen) < min_aln:
                continue
            key_qm = qmap.get(query)
            if key_qm is None:
                continue
            rank, clade, cluster = key_qm
            ri = RANK_POS.get(rank)
            if ri is None:
                continue
            try:
                lineage = idx_lineage.get(genome_idx(target))
            except (ValueError, IndexError):
                continue
            if lineage is None:
                continue
            target_clade = lineage[ri]
            if target_clade == clade or target_clade == "NA":
                continue  # in-clade hit (includes the marker's own copy)
            key = (rank, clade, cluster)
            ident = norm_identity(pident)
            n, max_id, worst = offenders.get(key, (0, 0.0, ""))
            if ident > max_id:
                max_id, worst = ident, target_clade
            offenders[key] = (n + 1, max_id, worst)
    return offenders


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True, help="mmseqs m8: query target pident alnlen")
    ap.add_argument("--idmap", required=True, help="m<N> <TAB> original query header")
    ap.add_argument("--markers", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--min_id", type=float, required=True)
    ap.add_argument("--min_aln", type=int, required=True)
    ap.add_argument("--out_markers", required=True)
    ap.add_argument("--out_fasta", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    idx_lineage = load_manifest(args.manifest)
    qmap = load_idmap(args.idmap)
    offenders = scan_hits(args.hits, qmap, idx_lineage, args.min_id, args.min_aln)

    kept = 0
    dropped = 0
    with open(args.markers) as fh, \
         open(args.out_markers, "w") as out, \
         open(args.report, "w") as rep:
        header = fh.readline().rstrip("\n")
        out.write(header + "\n")
        rep.write("rank\tclade\tcluster\tpass\tn_offtarget\tmax_offtarget_id\tworst_offtarget_clade\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rank, clade, cluster = f[0], f[1], f[2]
            key = (rank, clade, cluster)
            hit = offenders.get(key)
            if hit is None:
                out.write(line if line.endswith("\n") else line + "\n")
                rep.write(f"{rank}\t{clade}\t{cluster}\t1\t0\t0.0000\t\n")
                kept += 1
            else:
                n, max_id, worst = hit
                rep.write(f"{rank}\t{clade}\t{cluster}\t0\t{n}\t{max_id:.4f}\t{worst}\n")
                dropped += 1

    # Filter the FASTA to markers that survived (keyed by rank|clade|cluster).
    survivors = set()
    with open(args.markers) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[0], f[1], f[2])
            if key not in offenders:
                survivors.add(key)

    written = 0
    with open(args.fasta) as fh, open(args.out_fasta, "w") as out:
        emit = False
        for line in fh:
            if line.startswith(">"):
                parts = line[1:].rstrip("\n").split("|")
                key = (parts[0], parts[1], parts[2]) if len(parts) >= 3 else None
                emit = key in survivors
                if emit:
                    out.write(line)
                    written += 1
            elif emit:
                out.write(line)

    print(f"specificity_guard: kept {kept} markers, dropped {dropped} cross-mapping; "
          f"{written} sequences retained")


if __name__ == "__main__":
    main()
