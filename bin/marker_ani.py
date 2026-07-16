#!/usr/bin/env python3
"""
Per-species WITHIN-marker nucleotide distance, for the report's boxplot panels.

For each species we take its top-N markers by SCORE score and, for each marker,
compute the pairwise nucleotide distance among that marker's copies ACROSS THE
SPECIES' GENOMES -- i.e. how conserved the marker gene is within the species. A
low box = the marker is near-identical in every genome (a good, mappable marker);
a tall box / high whisker = it is polymorphic within the species, so reads from
divergent strains may not map.

The marker payload FASTA carries only one rep copy per marker, so the copies must
be gathered from the full per-genome CDS set:
    clusters_species.tsv   rep <TAB> member          -- marker cluster -> its members
    manifest.tsv           idx -> species            -- keep members in THIS species
    all_cds.ffn            >g<idx>_...  nucleotide    -- the member sequences

Distance is the mash estimate from k-mer Jaccard (same estimator as diversity.py,
expressed as distance = 1 - identity). Up to --max_copies genomes per marker are
used (subsampled in file order) to bound cost on large species.

Stages (each capped at --cap markers, x-axis ordered by score):
    specific-200   -- top-N species markers by score (pre cross-map guard)
    post-crossmap  -- those surviving the guard cleanly (pass, not mask-rescued)
    post-recovery  -- final set: clean survivors + mask-rescued

Output JSON embedded by build_report.py:
    {"kmer":..,"cap":..,"max_copies":..,"ymax":..,"stages":[...],
     "species": {sp: {stage: [{"c":cluster,"lo","q1","med","q3","hi","n":pairs}, ...]}}}
"""
import argparse
import json
import math
from itertools import combinations

STAGES = ["specific-200", "post-crossmap", "post-recovery"]


def strip_tag(c):
    return c[2:] if c[:2] in ("S:", "L:") else c


def genome_idx(seq_id):
    # g<idx>_<original token>
    return int(seq_id[1:].split("_", 1)[0])


def kmer_set(seq, k):
    seq = seq.upper()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else set()


def mash_dist(a, b, k):
    """Mash distance from k-mer Jaccard (0 = identical)."""
    if not a or not b:
        return 1.0
    j = len(a & b) / len(a | b)
    if j <= 0:
        return 1.0
    if j >= 1:
        return 0.0
    return max(0.0, -1.0 / k * math.log(2 * j / (1 + j)))


def load_manifest(path):
    idx_sp = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            idx_sp[int(f[ci["idx"]])] = f[ci["species"]]
    return idx_sp


def load_scored(path, rank):
    """species -> [(cluster_tagged, rep, score)] for the requested rank."""
    per_sp = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != rank:
                continue
            try:
                score = float(f[7]) if len(f) > 7 and f[7] != "" else 0.0
            except ValueError:
                score = 0.0
            per_sp.setdefault(f[1], []).append((f[2], strip_tag(f[2]), score))
    return per_sp


def load_specificity(path, rank):
    """(clade, cluster_tagged) -> (passed, recovered). Empty => no guard ran."""
    verdicts = {}
    if not path or path.endswith("NO_FILE"):
        return verdicts
    try:
        fh = open(path)
    except OSError:
        return verdicts
    with fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        rec_i = col.get("recovered")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] != rank:
                continue
            passed = f[3] == "1"
            recovered = rec_i is not None and rec_i < len(f) and f[rec_i] == "1"
            verdicts[(f[1], f[2])] = (passed, recovered)
    return verdicts


def quantile(sorted_vals, q):
    n = len(sorted_vals)
    if n == 1:
        return sorted_vals[0]
    pos = q * (n - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, n - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def box(vals):
    s = sorted(vals)
    return {
        "lo": round(s[0], 4), "q1": round(quantile(s, 0.25), 4),
        "med": round(quantile(s, 0.5), 4), "q3": round(quantile(s, 0.75), 4),
        "hi": round(s[-1], 4), "n": len(s),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scored", required=True, help="markers table (score col)")
    ap.add_argument("--specificity", help="specificity_report.tsv or NO_FILE")
    ap.add_argument("--clusters_species", required=True,
                    help="species clusters TSV: rep <TAB> member")
    ap.add_argument("--manifest", required=True, help="manifest.tsv (idx->species)")
    ap.add_argument("--all_cds", required=True,
                    help="all-genome CDS FASTA (nucleotide, >g<idx>_... headers)")
    ap.add_argument("--rank", default="species")
    ap.add_argument("--cap", type=int, default=200)
    ap.add_argument("--max_copies", type=int, default=60,
                    help="max genomes sampled per marker (bounds cost on big species)")
    ap.add_argument("--kmer", type=int, default=21)
    ap.add_argument("--ymax", type=float, default=0.2, help="boxplot y-axis max (distance)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    idx_sp = load_manifest(args.manifest)
    scored = load_scored(args.scored, args.rank)

    # Per species: top-cap markers by score. rep -> (species, cluster_tagged, order).
    rep_info = {}
    sp_order = {}   # species -> [rep,...] in score order (x-axis order)
    for sp, markers in scored.items():
        ranked = sorted(markers, key=lambda t: -t[2])[: args.cap]
        sp_order[sp] = [rep for _c, rep, _s in ranked]
        for ctag, rep, _s in ranked:
            rep_info[rep] = (sp, ctag)

    # Stream clusters_species: gather up to max_copies within-species members per
    # target marker. rep -> [member CDS ids]; also global set of needed CDS ids.
    marker_members = {}
    need = set()
    with open(args.clusters_species) as fh:
        for line in fh:
            i = line.find("\t")
            if i < 0:
                continue
            rep = line[:i]
            info = rep_info.get(rep)
            if info is None:
                continue
            member = line[i + 1:].rstrip("\n")
            lst = marker_members.setdefault(rep, [])
            if len(lst) >= args.max_copies:
                continue
            try:
                if idx_sp.get(genome_idx(member)) != info[0]:
                    continue  # member is not in this marker's species
            except (ValueError, IndexError):
                continue
            lst.append(member)
            need.add(member)

    # Stream all_cds once, pulling only the needed member sequences.
    seqs = {}
    cur, buf, keep = None, [], False
    with open(args.all_cds) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                keep = cur in need
                buf = []
            elif keep:
                buf.append(line.strip())
        if keep and cur is not None:
            seqs[cur] = "".join(buf)

    # Per marker: pairwise within-species distance among its copies.
    marker_box = {}   # rep -> box dict (distances) or None if < 2 copies
    for rep, members in marker_members.items():
        copies = [seqs[m] for m in members if m in seqs]
        if len(copies) < 2:
            continue
        ksets = [kmer_set(s, args.kmer) for s in copies]
        dists = [mash_dist(a, b, args.kmer) for a, b in combinations(ksets, 2)]
        if dists:
            marker_box[rep] = box(dists)

    verdicts = load_specificity(args.specificity, args.rank)
    have_qc = len(verdicts) > 0

    out = {"kmer": args.kmer, "cap": args.cap, "max_copies": args.max_copies,
           "ymax": args.ymax, "stages": STAGES, "species": {}}
    for sp, reps in sp_order.items():
        # only markers that actually have a within-species distance box.
        stage1 = [r for r in reps if r in marker_box]
        if not stage1:
            continue

        def keep_stage(pred):
            res = []
            for r in stage1:
                ctag = rep_info[r][1]
                v = verdicts.get((sp, ctag), (True, False))
                if pred(v):
                    res.append(r)
            return res

        if have_qc:
            s2 = keep_stage(lambda v: v[0] and not v[1])
            s3 = keep_stage(lambda v: v[0])
        else:
            s2 = s3 = list(stage1)
        sets = {"specific-200": stage1, "post-crossmap": s2, "post-recovery": s3}

        sp_out = {}
        for stage, members in sets.items():
            boxes = []
            for r in members:  # x-axis order = score order
                b = dict(marker_box[r])
                b["c"] = rep_info[r][1]
                boxes.append(b)
            sp_out[stage] = boxes
        out["species"][sp] = sp_out

    with open(args.out, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(f"marker_ani: {len(out['species'])} species, cap {args.cap}, "
          f"max_copies {args.max_copies}, k={args.kmer}, ymax {args.ymax}")


if __name__ == "__main__":
    main()
