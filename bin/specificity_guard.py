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
import os
import sys

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


# --------------------------------------------------- mask-recovery (k-mer) ----
# A cross-mapping marker can still be usable if the cross-mapping is confined to
# part of the gene: mask the shared stretches (k-mer matches to the off-target
# CDS) and, if a long enough clean window remains and not too much is masked, the
# marker is RECOVERED (kept whole -- the mask is only the recovery criterion).

_COMP = str.maketrans("ACGTacgt", "TGCATGCA")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def canonical_kmers(seq, k):
    seq = seq.upper()
    out = set()
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        rc = revcomp(km)
        out.add(km if km <= rc else rc)
    return out


def covered_positions(query, targets, k, merge_gap):
    """Merged intervals of query positions whose k-mer occurs in any target."""
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
    intervals, i = [], 0
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


def longest_clean(intervals, n):
    """Longest run of un-masked (clean) positions."""
    best, prev = 0, 0
    for s, e in intervals:
        best = max(best, s - prev)
        prev = e
    return max(best, n - prev)


def read_fasta_subset(path, wanted_ids):
    """Stream a FASTA once; {first_token_id: seq} for wanted ids only. Target
    headers are 'g<idx>_...' (no spaces), so the first whitespace token is the id."""
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


def scan_hits(path, qmap, idx_lineage, min_id, min_aln, collect_targets=False):
    """Return {(rank,clade,cluster): (n_offtarget, max_id, worst_clade, {clades})}
    and, when collect_targets, {key: {off_target_gene_id}} for mask recovery.

    The clade set is EVERY off-target clade the marker's nucleotide sequence hits
    (at the marker's rank), not just the worst -- so a downstream merge probe can
    tell whether pooling some of those clades into the marker's clade would make it
    unique (all off-targets absorbed) or not (an off-target remains outside)."""
    offenders = {}
    off_by_clade = {}      # key -> {off_target_clade: {target_gene_id}}
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
            n, max_id, worst, clades = offenders.get(key, (0, 0.0, "", set()))
            if ident > max_id:
                max_id, worst = ident, target_clade
            clades.add(target_clade)
            offenders[key] = (n + 1, max_id, worst, clades)
            if collect_targets:
                off_by_clade.setdefault(key, {}).setdefault(
                    target_clade, set()).add(target)
    return offenders, off_by_clade


def marker_rep_seqs(fasta_path, wanted_keys):
    """First CDS sequence per (rank,clade,cluster) among wanted_keys. Marker
    headers are 'rank|clade|cluster|genome_id' (clade may contain spaces, so key
    off the '|' split, not whitespace)."""
    seqs, cur, buf = {}, None, []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None and cur not in seqs:
                    seqs[cur] = "".join(buf)
                p = line[1:].rstrip("\n").split("|")
                k = (p[0], p[1], p[2]) if len(p) >= 3 else None
                cur = k if (k in wanted_keys and k not in seqs) else None
                buf = []
            elif cur is not None:
                buf.append(line.strip())
        if cur is not None and cur not in seqs:
            seqs[cur] = "".join(buf)
    return seqs


def union_len(intervals):
    return sum(e - s for s, e in intervals)


def merge_intervals(intervals):
    """Sort + merge overlapping [s,e) intervals."""
    if not intervals:
        return []
    iv = sorted(intervals)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def recover_markers(offenders, off_by_clade, args):
    """Which cross-mapping markers are recoverable by masking: a clean (un-cross-
    mapping) window of >= recovery_min_clean bp remains AND <= recovery_max_masked
    fraction of the gene is masked.

    Returns
      verdict: {key: (clean_bp, masked_frac, bool)}   -- recovery vs ALL off-targets
      by_clade: {key: (gene_len, {off_clade: [(s,e),...]})}  -- per-off-target-clade
                masked intervals, so a merge probe can re-mask against only the
                clades that stay OUTSIDE a merged clade (merging absorbs the rest)."""
    all_targets = set()
    for cm in off_by_clade.values():
        for gs in cm.values():
            all_targets |= gs
    sys.stderr.write(f"[guard] mask recovery: {len(offenders)} dropped markers, "
                     f"{len(all_targets)} off-target genes; loading targets from "
                     f"{os.path.basename(args.target)}...\n")
    target_seqs = read_fasta_subset(args.target, all_targets)
    marker_seqs = marker_rep_seqs(args.fasta, set(offenders))

    verdict, by_clade = {}, {}
    for key in offenders:
        qseq = marker_seqs.get(key, "")
        n = len(qseq)
        cmap = off_by_clade.get(key, {})
        per_clade = {}
        for oc, gene_ids in cmap.items():
            tgts = [target_seqs[t] for t in gene_ids if t in target_seqs]
            iv = covered_positions(qseq, tgts, args.k, args.merge_gap) if n else []
            if iv:
                per_clade[oc] = iv
        by_clade[key] = (n, per_clade)
        # Recovery vs ALL off-target clades = union of the per-clade intervals.
        allmask = merge_intervals([iv for ivs in per_clade.values() for iv in ivs])
        clean = longest_clean(allmask, n) if n else 0
        frac = (union_len(allmask) / n) if n else 1.0
        rec = bool(n) and clean >= args.recovery_min_clean \
            and frac <= args.recovery_max_masked_frac
        verdict[key] = (clean, frac, rec)
    return verdict, by_clade


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
    ap.add_argument("--out_mask_intervals",
                    help="per-off-target-clade masked intervals (for re-masking "
                    "after a hypothetical species merge)")
    # Mask recovery: keep a cross-mapping marker if masking the shared stretches
    # still leaves a long clean window. --target is the off-target CDS universe
    # (all_cds.ffn). Omit --target to disable recovery (legacy binary guard).
    ap.add_argument("--target", help="off-target CDS (all_cds.ffn) for mask recovery")
    ap.add_argument("--recovery_min_clean", type=int, default=300,
                    help="min clean (un-cross-mapping) window bp to recover")
    ap.add_argument("--recovery_max_masked_frac", type=float, default=0.5,
                    help="max fraction of the gene that may be masked")
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--merge_gap", type=int, default=50)
    args = ap.parse_args()

    idx_lineage = load_manifest(args.manifest)
    qmap = load_idmap(args.idmap)
    do_recover = bool(args.target)
    offenders, off_by_clade = scan_hits(args.hits, qmap, idx_lineage,
                                        args.min_id, args.min_aln,
                                        collect_targets=do_recover)
    if do_recover:
        recovered, by_clade = recover_markers(offenders, off_by_clade, args)
    else:
        recovered, by_clade = {}, {}

    # Per-off-target-clade masked intervals, so a merge probe can re-mask a marker
    # against only the clades that remain OUTSIDE a merged clade.
    if args.out_mask_intervals:
        with open(args.out_mask_intervals, "w") as mi:
            mi.write("rank\tclade\tcluster\tgene_len\tclade_masks\n")
            for key, (glen, per_clade) in by_clade.items():
                enc = "|".join(
                    f"{oc}:" + ",".join(f"{s}-{e}" for s, e in ivs)
                    for oc, ivs in per_clade.items() if ivs)
                mi.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{glen}\t{enc}\n")

    kept = dropped = rescued = 0
    survivors = set()
    with open(args.markers) as fh, \
         open(args.out_markers, "w") as out, \
         open(args.report, "w") as rep:
        header = fh.readline().rstrip("\n")
        out.write(header + "\n")
        rep.write("rank\tclade\tcluster\tpass\tn_offtarget\tmax_offtarget_id\t"
                  "worst_offtarget_clade\tofftarget_clades\trecovered\t"
                  "clean_bp\tmasked_frac\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rank, clade, cluster = f[0], f[1], f[2]
            key = (rank, clade, cluster)
            hit = offenders.get(key)
            if hit is None:
                # No off-target at all: clean marker.
                out.write(line if line.endswith("\n") else line + "\n")
                rep.write(f"{rank}\t{clade}\t{cluster}\t1\t0\t0.0000\t\t\t0\t\t\n")
                survivors.add(key)
                kept += 1
                continue
            n, max_id, worst, clades = hit
            off = ";".join(sorted(clades))
            clean, frac, rec = recovered.get(key, (0, 1.0, False))
            # pass reflects the FINAL kept status (clean OR recovered) so every
            # downstream reader treats recovered markers as final markers.
            passed = 1 if rec else 0
            rep.write(f"{rank}\t{clade}\t{cluster}\t{passed}\t{n}\t{max_id:.4f}\t"
                      f"{worst}\t{off}\t{1 if rec else 0}\t{clean}\t{frac:.4f}\n")
            if rec:
                out.write(line if line.endswith("\n") else line + "\n")
                survivors.add(key)
                rescued += 1
            else:
                dropped += 1

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

    print(f"specificity_guard: kept {kept} clean + {rescued} mask-recovered markers, "
          f"dropped {dropped} cross-mapping; {written} sequences retained")


if __name__ == "__main__":
    main()
