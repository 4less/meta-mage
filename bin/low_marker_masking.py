#!/usr/bin/env python3
"""
Batch crossmap-masking diagnostics for low-marker species (pipeline driver).

Same analysis as assess_crossmap_masking.py, but driven from the marker tables:
it finds every species that ended with FEWER than --threshold final markers and,
for each, asks whether the markers the cross-map guard dropped are wholly
cross-mapping (masking can't help -> real duplicate / merge case) or only
partially (mask the offending slice and keep a clade-unique marker).

All flagged species are processed in ONE pass over each big input (crossmap.m8,
the target CDS, the marker FASTA), so cost is independent of how many species are
flagged. --target must be the SAME corpus CROSSMAP searched (all_cds.ffn): the
off-target ids come from crossmap.m8, so a narrower target (e.g. reps.ffn) cannot
resolve hits to non-rep genomes and silently writes them off as
'target_unavailable'. Only the off-target hits are held in memory, not the corpus.

Outputs (under --outdir):
    low_marker_species.tsv  -- species, genus, genomes, final_markers, flagged
                               (drives the ANI-gap step and the report)
    masking_summary.tsv     -- per flagged species: dropped, whole_gene, rescuable,
                               ambiguous, target_unavailable, pct_rescuable
    masking_markers.tsv     -- per dropped marker: spans, mask ranges, verdict
    masking_coverage.txt    -- ASCII cross-map track per marker
"""
import argparse
import os
import sys
from collections import defaultdict

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_POS = {r: i for i, r in enumerate(RANKS)}
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def load_manifest(path):
    """idx->lineage tuple; species->genus; species->genome count."""
    idx_lineage = {}
    sp_genus = {}
    sp_size = defaultdict(int)
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            idx_lineage[int(f[col["idx"]])] = tuple(f[col[r]] for r in RANKS)
            sp = f[col["species"]]
            if sp != "NA":
                sp_size[sp] += 1
                sp_genus[sp] = f[col["genus"]]
    return idx_lineage, sp_genus, sp_size


def load_idmap(path):
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


def final_counts(markers_specific):
    """species -> # final (post-QC) markers, from markers.specific.tsv."""
    cnt = defaultdict(int)
    with open(markers_specific) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0] == "species":
                cnt[f[1]] += 1
    return cnt


def dropped_by_species(report):
    """species -> set(cluster) dropped by the guard (rank=species, pass=0)."""
    out = defaultdict(set)
    with open(report) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[0] == "species" and f[3] == "0":
                out[f[1]].add(f[2])
    return out


def genome_idx(seq_id):
    return int(seq_id[1:].split("_", 1)[0])


def norm_id(v):
    v = float(v)
    return v / 100.0 if v > 1.0 else v


def read_fasta_subset(path, wanted_ids):
    """Stream a FASTA once; return {first_token_id: seq} for wanted ids only."""
    want = set(wanted_ids)
    seqs, cur, buf, emit = {}, None, [], False
    if not want:
        return seqs
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if emit and cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
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
    best, prev = 0, 0
    for s, e in intervals:
        best = max(best, s - prev)
        prev = e
    return max(best, n - prev)


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
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--markers", required=True, help="markers.specific.tsv (final)")
    ap.add_argument("--report", required=True, help="specificity_report.tsv")
    ap.add_argument("--idmap", required=True, help="query.map")
    ap.add_argument("--crossmap", required=True, help="crossmap.m8")
    ap.add_argument("--marker_fasta", required=True, help="markers.nuc.fasta")
    ap.add_argument("--target", required=True, help="target CDS (reps.ffn)")
    ap.add_argument("--threshold", type=int, default=50)
    ap.add_argument("--min_id", type=float, default=0.95)
    ap.add_argument("--min_aln", type=int, default=100)
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--merge_gap", type=int, default=50)
    ap.add_argument("--cover_frac", type=float, default=0.90)
    ap.add_argument("--min_survivor", type=int, default=200)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    idx_lineage, sp_genus, sp_size = load_manifest(args.manifest)
    q2key, key2q = load_idmap(args.idmap)
    fcount = final_counts(args.markers)
    drops = dropped_by_species(args.report)

    # Flag species below threshold (0 disables everything).
    flagged = sorted(sp for sp in sp_size
                     if args.threshold and fcount.get(sp, 0) < args.threshold)
    os.makedirs(args.outdir, exist_ok=True)

    with open(os.path.join(args.outdir, "low_marker_species.tsv"), "w") as fh:
        fh.write("species\tgenus\tgenomes\tfinal_markers\tflagged\n")
        for sp in sorted(sp_size):
            fl = "yes" if sp in set(flagged) else "no"
            fh.write(f"{sp}\t{sp_genus.get(sp,'NA')}\t{sp_size[sp]}\t"
                     f"{fcount.get(sp,0)}\t{fl}\n")

    # m-id -> (species, cluster) for flagged species' dropped markers.
    wanted_q = {}
    for sp in flagged:
        for cl in drops.get(sp, ()):
            q = key2q.get(("species", sp, cl))
            if q:
                wanted_q[q] = (sp, cl)
    if not wanted_q:
        sys.stderr.write("[lowmark] no dropped markers among flagged species; "
                         "writing empty masking outputs.\n")
        open(os.path.join(args.outdir, "masking_summary.tsv"), "w").write(
            "species\tgenomes\tfinal_markers\tdropped\twhole_gene\trescuable\t"
            "ambiguous\ttarget_unavailable\tpct_rescuable\n")
        open(os.path.join(args.outdir, "masking_markers.tsv"), "w").write(
            "species\tcluster\tgene_len\tn_offtarget_genes\tcovered_bp\t"
            "covered_frac\tlongest_safe_window\tmask_ranges_1based\tverdict\n")
        open(os.path.join(args.outdir, "masking_coverage.txt"), "w").write("")
        return

    # One pass over crossmap.m8: off-target target genes per marker.
    off = defaultdict(set)
    with open(args.crossmap) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            q, t, pid, aln = f[0], f[1], f[2], f[3]
            if q not in wanted_q:
                continue
            if norm_id(pid) < args.min_id or int(aln) < args.min_aln:
                continue
            rank, clade, _cl = q2key[q]
            ri = RANK_POS[rank]
            try:
                lineage = idx_lineage.get(genome_idx(t))
            except (ValueError, IndexError):
                continue
            if lineage is None:
                continue
            tclade = lineage[ri]
            if tclade == clade or tclade == "NA":
                continue
            off[q].add(t)

    all_targets = set().union(*off.values()) if off else set()
    sys.stderr.write(f"[lowmark] {len(flagged)} flagged species, {len(wanted_q)} "
                     f"dropped markers, {len(all_targets)} off-target genes; "
                     f"scanning {os.path.basename(args.target)} ...\n")
    target_seqs = read_fasta_subset(args.target, all_targets)

    # One pass over the marker FASTA: query seq per (species, cluster) we need.
    need_clusters = {(sp, cl) for sp, cl in wanted_q.values()}
    marker_seqs = {}
    keep, buf = None, []
    with open(args.marker_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep is not None:
                    marker_seqs[keep] = "".join(buf)
                p = line[1:].rstrip("\n").split("|")
                keep = (p[1], p[2]) if (len(p) >= 3 and p[0] == "species"
                                        and (p[1], p[2]) in need_clusters) else None
                buf = []
            elif keep is not None:
                buf.append(line.strip())
        if keep is not None:
            marker_seqs[keep] = "".join(buf)

    # Score every dropped marker.
    per_species = defaultdict(lambda: defaultdict(int))
    marker_rows = []
    tracks = []
    for q, (sp, cl) in sorted(wanted_q.items()):
        qseq = marker_seqs.get((sp, cl), "")
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
            verdict = "target_unavailable"
        elif frac >= args.cover_frac:
            verdict = "whole_gene"
        elif safe_win >= args.min_survivor:
            verdict = "rescuable"
        else:
            verdict = "ambiguous"
        per_species[sp][verdict] += 1
        per_species[sp]["dropped"] += 1
        mask_ranges = ";".join(f"{s+1}-{e}" for s, e in intervals) or "-"
        marker_rows.append((sp, cl, n, len(tgt_ids), covered_bp, round(frac, 3),
                            safe_win, mask_ranges, verdict))
        if n:
            tracks.append(f">{sp} | {cl}  len={n}bp  covered={frac:.0%}  "
                          f"safe_window={safe_win}bp  [{verdict}]\n"
                          f"  {ascii_track(intervals, n)}\n")

    marker_rows.sort(key=lambda r: (r[0], r[5], -r[6]))
    with open(os.path.join(args.outdir, "masking_markers.tsv"), "w") as fh:
        fh.write("species\tcluster\tgene_len\tn_offtarget_genes\tcovered_bp\t"
                 "covered_frac\tlongest_safe_window\tmask_ranges_1based\tverdict\n")
        for r in marker_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    with open(os.path.join(args.outdir, "masking_summary.tsv"), "w") as fh:
        fh.write("species\tgenomes\tfinal_markers\tdropped\twhole_gene\trescuable\t"
                 "ambiguous\ttarget_unavailable\tpct_rescuable\n")
        for sp in flagged:
            s = per_species.get(sp)
            if not s:
                continue
            d = s.get("dropped", 0)
            resc = s.get("rescuable", 0)
            pct = f"{100 * resc / d:.0f}" if d else "0"
            fh.write(f"{sp}\t{sp_size[sp]}\t{fcount.get(sp,0)}\t{d}\t"
                     f"{s.get('whole_gene',0)}\t{resc}\t{s.get('ambiguous',0)}\t"
                     f"{s.get('target_unavailable',0)}\t{pct}\n")

    with open(os.path.join(args.outdir, "masking_coverage.txt"), "w") as fh:
        fh.write(f"# k={args.k} merge_gap={args.merge_gap} min_id={args.min_id} "
                 f"min_aln={args.min_aln} min_survivor={args.min_survivor} "
                 f"(# = cross-maps off-target)\n\n")
        fh.writelines(tracks)

    tot = len(marker_rows)
    resc = sum(1 for r in marker_rows if r[8] == "rescuable")
    sys.stderr.write(f"[lowmark] {tot} dropped markers assessed; {resc} rescuable "
                     f"by masking. -> {args.outdir}\n")


if __name__ == "__main__":
    main()
