#!/usr/bin/env python3
"""
Per-species pairwise marker ANI, for the report's boxplot panels.

For each species we take its top-N markers by SCORE score and, at three filtering
stages, compute the pairwise nucleotide identity BETWEEN markers (mash-style k-mer
identity, same estimator as diversity.py). Each marker's box summarises its
identity to every other marker in the stage's set -- a high box means the marker
is redundant with the rest of the set, a low box means it is distinct.

Stages (each capped at --cap markers, x-axis ordered by score):
    specific-200   -- top-N species markers by score (pre cross-map guard)
    post-crossmap  -- those surviving the guard cleanly (pass, not mask-rescued)
    post-recovery  -- final set: clean survivors + mask-rescued

Output is compact JSON embedded by build_report.py:
    {"kmer":..., "cap":..., "stages":[...],
     "species": {sp: {stage: [{"c":cluster,"lo","q1","med","q3","hi","n":pairs}, ...]}}}
"""
import argparse
import json
import math
from itertools import combinations

STAGES = ["specific-200", "post-crossmap", "post-recovery"]


def kmer_set(seq, k):
    seq = seq.upper()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else set()


def mash_identity(a, b, k):
    if not a or not b:
        return 0.0
    j = len(a & b) / len(a | b)
    if j <= 0:
        return 0.0
    if j >= 1:
        return 1.0
    d = -1.0 / k * math.log(2 * j / (1 + j))
    return max(0.0, 1.0 - d)


def read_rep_seqs(fasta):
    """(rank,clade,cluster) -> first CDS sequence (one rep copy per marker)."""
    seqs, key, buf = {}, None, []

    def flush():
        if key is not None and key not in seqs:
            seqs[key] = "".join(buf)

    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                p = line[1:].strip().split("|")
                key = (p[0], p[1], p[2]) if len(p) >= 3 else None
                buf = []
            else:
                buf.append(line.strip())
        flush()
    return seqs


def load_scored(path, rank):
    """species clade -> [(cluster, score)] for the requested rank."""
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
            per_sp.setdefault(f[1], []).append((f[2], score))
    return per_sp


def load_specificity(path, rank):
    """(clade,cluster) -> (passed, recovered). Empty dict => no guard ran."""
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
    """Linear-interpolation percentile on a pre-sorted list."""
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
    ap.add_argument("--fasta", required=True, help="marker CDS FASTA (nucleotide)")
    ap.add_argument("--scored", required=True, help="markers.specific.tsv (score col)")
    ap.add_argument("--specificity", help="specificity_report.tsv or NO_FILE")
    ap.add_argument("--rank", default="species")
    ap.add_argument("--cap", type=int, default=200)
    ap.add_argument("--kmer", type=int, default=21)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rep = read_rep_seqs(args.fasta)
    scored = load_scored(args.scored, args.rank)
    verdicts = load_specificity(args.specificity, args.rank)
    have_qc = len(verdicts) > 0

    out = {"kmer": args.kmer, "cap": args.cap, "stages": STAGES, "species": {}}
    for sp, markers in scored.items():
        # stage1: top-cap by score, only markers with a sequence.
        ranked = sorted(markers, key=lambda cs: -cs[1])
        stage1 = [c for c, _ in ranked if (args.rank, sp, c) in rep][: args.cap]
        if len(stage1) < 2:
            continue
        ksets = {c: kmer_set(rep[(args.rank, sp, c)], args.kmer) for c in stage1}
        # full pairwise identity among stage1 markers (reused by all stages).
        idm = {}
        for a, b in combinations(stage1, 2):
            v = mash_identity(ksets[a], ksets[b], args.kmer)
            idm[(a, b)] = v
            idm[(b, a)] = v

        def stage_set(pred):
            return [c for c in stage1 if pred(c)]

        if have_qc:
            s2 = stage_set(lambda c: verdicts.get((sp, c), (True, False))[0]
                           and not verdicts.get((sp, c), (True, False))[1])
            s3 = stage_set(lambda c: verdicts.get((sp, c), (True, False))[0])
        else:
            s2 = s3 = list(stage1)
        sets = {"specific-200": stage1, "post-crossmap": s2, "post-recovery": s3}

        sp_out = {}
        for stage, members in sets.items():
            mset = set(members)
            boxes = []
            for c in members:  # x-axis order = score order (stage1 order preserved)
                dists = [idm[(c, o)] for o in members if o != c]
                if not dists:
                    continue
                b = box(dists)
                b["c"] = c
                boxes.append(b)
            sp_out[stage] = boxes
        out["species"][sp] = sp_out

    with open(args.out, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(f"marker_ani: {len(out['species'])} species with >=2 markers, "
          f"cap {args.cap}, k={args.kmer}")


if __name__ == "__main__":
    main()
