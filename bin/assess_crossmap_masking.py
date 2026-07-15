#!/usr/bin/env python3
"""
Does cross-mapping cover a whole marker, or could masking part of it rescue it?

The specificity guard drops a marker outright the moment its nucleotide sequence
maps to ONE off-target genome over a read-length window. But if that cross-mapping
region is only a slice of the gene, the rest may still be clade-unique -- masking
the offending slice (MetaPhlAn stores markers with masked ranges) would keep a
usable, specific marker instead of discarding it.

For one species this tool takes every marker the guard dropped and reconstructs,
per marker, WHERE on the gene the off-target hits land:

  1. find the species' dropped markers in specificity_report.tsv (pass == 0);
  2. read crossmap.m8 for those markers' off-target target genes (same off-target
     test and id/aln thresholds the guard used);
  3. pull the marker CDS (query) and those off-target CDS (target) sequences;
  4. mark the query positions homologous to ANY off-target gene via shared
     canonical k-mers (strand-agnostic; small gaps between anchors are filled),
     giving the cross-mapped span vs the safe (clade-unique) span;
  5. classify each marker:
       whole_gene  -- >= --cover_frac of the gene cross-maps; masking won't help
                      (this is a real duplicate region / merge-or-drop case);
       rescuable   -- the safe complement has a contiguous window >= --min_survivor;
                      mask the cross-mapped range(s) and keep the marker;
       ambiguous   -- cross-mapped but the largest safe window is too short to
                      trust for read mapping.

Outputs (under --outdir):
    <species>.masking.tsv   -- one row per dropped marker with spans + verdict +
                               suggested mask interval(s) on the query (1-based)
    <species>.coverage.txt  -- an ASCII cross-map track per marker for eyeballing

k-mer overlap is an alignment-free proxy for the aligned region; it needs no
mmseqs/skani, only the run's own files. The only heavy input is emit/all_cds.ffn
(all genomes' CDS), streamed once to pull the handful of target genes we need.
"""
import argparse
import os
import subprocess
import sys
from collections import defaultdict

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_POS = {r: i for i, r in enumerate(RANKS)}
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def load_manifest(path):
    """idx -> lineage tuple (one clade per rank)."""
    idx_lineage = {}
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            idx_lineage[int(f[col["idx"]])] = tuple(f[col[r]] for r in RANKS)
    return idx_lineage


def load_idmap(path):
    """m<N> -> (rank, clade, cluster). Also build (rank,clade,cluster) -> m<N>."""
    q2key, key2q = {}, {}
    with open(path) as fh:
        for line in fh:
            qid, header = line.rstrip("\n").split("\t", 1)
            parts = header.split("|")
            if len(parts) >= 3:
                key = (parts[0], parts[1], parts[2])
                q2key[qid] = key
                key2q[key] = qid
    return q2key, key2q


def dropped_clusters(report, species):
    """clusters (with S:/L: tag) the guard dropped for this species (rank=species)."""
    out = set()
    with open(report) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            rank, clade, cluster, passed = f[0], f[1], f[2], f[3]
            if rank == "species" and clade == species and passed == "0":
                out.add(cluster)
    return out


def genome_idx(seq_id):
    return int(seq_id[1:].split("_", 1)[0])


def norm_id(v):
    v = float(v)
    return v / 100.0 if v > 1.0 else v


def collect_offtarget_targets(crossmap, wanted_q, q2key, idx_lineage, min_id, min_aln):
    """m-id -> set(target_gene_id) for off-target hits passing the guard's cuts."""
    hits = defaultdict(set)
    with open(crossmap) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            q, t, pid, aln = f[0], f[1], f[2], f[3]
            if q not in wanted_q:
                continue
            if norm_id(pid) < min_id or int(aln) < min_aln:
                continue
            key = q2key.get(q)
            if not key:
                continue
            rank, clade, _cl = key
            ri = RANK_POS[rank]
            try:
                lineage = idx_lineage.get(genome_idx(t))
            except (ValueError, IndexError):
                continue
            if lineage is None:
                continue
            tclade = lineage[ri]
            if tclade == clade or tclade == "NA":
                continue                       # in-clade: ignored, same as guard
            hits[q].add(t)
    return hits


def read_fasta_subset(path, wanted_ids, by_first_token=True):
    """Stream a (possibly huge) FASTA once; return {id: seq} for wanted ids only."""
    want = set(wanted_ids)
    seqs, cur, buf, emit = {}, None, [], False
    if not want:
        return seqs
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if emit and cur is not None:
                    seqs[cur] = "".join(buf)
                header = line[1:].rstrip("\n")
                cur = header.split()[0] if by_first_token else header
                emit = cur in want
                buf = []
            elif emit:
                buf.append(line.strip())
        if emit and cur is not None:
            seqs[cur] = "".join(buf)
    return seqs


def canonical_kmers(seq, k):
    seq = seq.upper()
    out = set()
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        rc = revcomp(km)
        out.add(km if km <= rc else rc)
    return out


def covered_positions(query, targets, k, merge_gap):
    """Return merged [start,end) intervals (0-based) of query covered by any target.

    A query position i is 'covered' if the k-mer starting at i (canonical) occurs
    in any target. Runs are then merged across gaps <= merge_gap so a homologous
    region with sparse anchors reads as one interval.
    """
    q = query.upper()
    n = len(q)
    if n < k or not targets:
        return []
    tk = set()
    for t in targets:
        tk |= canonical_kmers(t, k)
    covered = bytearray(n)
    for i in range(n - k + 1):
        km = q[i:i + k]
        rc = revcomp(km)
        if (km if km <= rc else rc) in tk:
            for j in range(i, i + k):
                covered[j] = 1
    # Merge into intervals with small-gap filling.
    intervals = []
    i = 0
    while i < n:
        if covered[i]:
            j = i
            while j < n and covered[j]:
                j += 1
            intervals.append([i, j])
            i = j
        else:
            i += 1
    if not intervals:
        return []
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        if s - merged[-1][1] <= merge_gap:
            merged[-1][1] = e
        else:
            merged.append([s, e])
    return merged


def longest_gap(intervals, n):
    """Longest contiguous uncovered window given covered intervals over [0,n)."""
    best, prev = 0, 0
    for s, e in intervals:
        best = max(best, s - prev)
        prev = e
    best = max(best, n - prev)
    return best


def ascii_track(intervals, n, width=60):
    cells = ["."] * width
    for s, e in intervals:
        a = int(s / n * width)
        b = max(a + 1, int(e / n * width))
        for i in range(a, min(b, width)):
            cells[i] = "#"
    return "".join(cells)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True, help='e.g. "s__Bacteroides caecimuris"')
    ap.add_argument("--base", required=True, help="nextflow outdir")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min_id", type=float, default=0.95,
                    help="off-target identity cut (match specificity_min_id)")
    ap.add_argument("--min_aln", type=int, default=100,
                    help="off-target min alignment bp (match specificity_min_aln)")
    ap.add_argument("--k", type=int, default=21, help="k-mer size (default 21)")
    ap.add_argument("--merge_gap", type=int, default=50,
                    help="fill gaps <= this between anchors (default 50 bp)")
    ap.add_argument("--cover_frac", type=float, default=0.90,
                    help=">= this covered fraction => whole gene cross-maps (0.90)")
    ap.add_argument("--min_survivor", type=int, default=200,
                    help="min contiguous safe window to call a marker rescuable (bp)")
    ap.add_argument("--target_cds", help="override target FASTA "
                    "(default <base>/emit/all_cds.ffn)")
    ap.add_argument("--markers_fasta", help="override marker CDS "
                    "(default <base>/markers/markers.nuc.fasta)")
    args = ap.parse_args()

    b = args.base
    manifest = os.path.join(b, "manifest", "manifest.tsv")
    report = os.path.join(b, "markers", "specificity_report.tsv")
    idmap = os.path.join(b, "markers", "query.map")
    crossmap = os.path.join(b, "markers", "crossmap.m8")
    markers_fa = args.markers_fasta or os.path.join(b, "markers", "markers.nuc.fasta")
    target_fa = args.target_cds or os.path.join(b, "emit", "all_cds.ffn")
    for pth in (manifest, report, idmap, crossmap, markers_fa, target_fa):
        if not os.path.exists(pth):
            sys.exit(f"error: missing expected input: {pth}")

    idx_lineage = load_manifest(manifest)
    q2key, key2q = load_idmap(idmap)
    drop = dropped_clusters(report, args.species)
    if not drop:
        sys.exit(f"{args.species}: no guard-dropped markers in {report} "
                 "(nothing to assess).")

    # m-ids of this species' dropped markers.
    wanted_q, q2cluster = {}, {}
    for cl in drop:
        q = key2q.get(("species", args.species, cl))
        if q:
            wanted_q[q] = cl
            q2cluster[q] = cl
    if not wanted_q:
        sys.exit(f"{args.species}: {len(drop)} dropped clusters but none found in "
                 f"query.map -- id-map/species mismatch?")
    sys.stderr.write(f"[mask] {args.species}: {len(wanted_q)} dropped markers\n")

    off = collect_offtarget_targets(crossmap, set(wanted_q), q2key, idx_lineage,
                                    args.min_id, args.min_aln)
    all_targets = set().union(*off.values()) if off else set()
    sys.stderr.write(f"[mask] {len(all_targets)} distinct off-target genes; "
                     f"scanning {os.path.basename(target_fa)} ...\n")

    # Marker CDS: headers are rank|clade|cluster|genome; key surviving rows by cluster.
    marker_seqs = {}   # cluster -> seq
    cur, buf, keep_cl = None, [], None
    with open(markers_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep_cl is not None:
                    marker_seqs[keep_cl] = "".join(buf)
                parts = line[1:].rstrip("\n").split("|")
                keep_cl = parts[2] if (len(parts) >= 3 and parts[0] == "species"
                                       and parts[1] == args.species
                                       and parts[2] in drop) else None
                buf = []
            elif keep_cl is not None:
                buf.append(line.strip())
        if keep_cl is not None:
            marker_seqs[keep_cl] = "".join(buf)

    target_seqs = read_fasta_subset(target_fa, all_targets)
    sys.stderr.write(f"[mask] pulled {len(target_seqs)}/{len(all_targets)} target "
                     f"seqs, {len(marker_seqs)} marker seqs\n")

    os.makedirs(args.outdir, exist_ok=True)
    safe_name = args.species.replace(" ", "_").replace("/", "_")
    tsv_path = os.path.join(args.outdir, f"{safe_name}.masking.tsv")
    trk_path = os.path.join(args.outdir, f"{safe_name}.coverage.txt")

    verdict_n = defaultdict(int)
    rows = []
    tracks = []
    for q, cl in sorted(wanted_q.items()):
        qseq = marker_seqs.get(cl, "")
        n = len(qseq)
        tgt_ids = off.get(q, set())
        tgts = [target_seqs[t] for t in tgt_ids if t in target_seqs]
        intervals = covered_positions(qseq, tgts, args.k, args.merge_gap) if n else []
        covered_bp = sum(e - s for s, e in intervals)
        frac = covered_bp / n if n else 0.0
        safe_win = longest_gap(intervals, n) if n else 0
        if n == 0:
            verdict = "no_seq"
        elif not tgts:
            # The marker cross-maps (it's in `off`) but none of its off-target
            # genes were in the target FASTA -- e.g. reps.ffn holds only rep
            # genomes and every hit here was to a non-rep genome. Can't localise
            # without a full all_cds.ffn pass; don't mistake this for "clean".
            verdict = "target_unavailable"
        elif frac >= args.cover_frac:
            verdict = "whole_gene"
        elif safe_win >= args.min_survivor:
            verdict = "rescuable"
        else:
            verdict = "ambiguous"
        verdict_n[verdict] += 1
        mask_ranges = ";".join(f"{s + 1}-{e}" for s, e in intervals) or "-"
        rows.append((cl, n, len(tgt_ids), covered_bp, round(frac, 3),
                     safe_win, mask_ranges, verdict))
        if n:
            tracks.append(f">{cl}  len={n}bp  covered={frac:.0%}  "
                          f"safe_window={safe_win}bp  [{verdict}]\n"
                          f"  {ascii_track(intervals, n)}  (# = cross-maps off-target)\n")

    rows.sort(key=lambda r: (r[4], -r[5]))  # most-covered / smallest-safe first
    with open(tsv_path, "w") as fh:
        fh.write("cluster\tgene_len\tn_offtarget_genes\tcovered_bp\tcovered_frac\t"
                 "longest_safe_window\tmask_ranges_1based\tverdict\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(trk_path, "w") as fh:
        fh.write(f"# {args.species}: cross-map coverage of guard-dropped markers\n")
        fh.write(f"# k={args.k} merge_gap={args.merge_gap} min_id={args.min_id} "
                 f"min_aln={args.min_aln} min_survivor={args.min_survivor}\n\n")
        fh.writelines(tracks)

    total = sum(verdict_n.values())
    sys.stderr.write(
        f"[mask] {args.species}: {total} dropped markers -> "
        + ", ".join(f"{k}={v}" for k, v in sorted(verdict_n.items())) + "\n")
    sys.stderr.write(f"[mask] wrote {tsv_path}\n[mask] wrote {trk_path}\n")
    n_resc = verdict_n.get("rescuable", 0)
    if total:
        sys.stderr.write(
            f"[mask] => {n_resc}/{total} ({100*n_resc/total:.0f}%) could be kept by "
            f"masking the cross-mapped range instead of dropping the whole marker.\n")


if __name__ == "__main__":
    main()
