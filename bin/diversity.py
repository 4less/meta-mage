#!/usr/bin/env python3
"""
Per-marker within-clade nucleotide divergence, from the species-rep copies.

Groups sequences by (rank, clade, cluster) from the header
    >{rank}|{clade}|{cluster}|{genome_id}
and computes pairwise Mash-style nucleotide identity among the copies of each
marker within its clade. Augments the scored marker table with:
    n_rep_copies  mean_within_id  min_within_id

Implementation is stdlib-only exact k-mer Jaccard (each marker has only a handful
of rep copies, so exact sets are cheap and accurate). At full scale swap in
`mash sketch -i` + `mash dist` per marker block; the columns are identical.
"""
import argparse
import math
from collections import defaultdict
from itertools import combinations


def read_groups(fasta):
    groups = defaultdict(list)  # (rank,clade,cluster) -> [seq,...]
    key = None
    buf = []

    def flush():
        if key is not None:
            groups[key[:3]].append(("".join(buf), key[3]))

    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                parts = line[1:].strip().split("|")
                # rank, clade, cluster, genome_id
                key = (parts[0], parts[1], parts[2], parts[3] if len(parts) > 3 else "")
                buf = []
            else:
                buf.append(line.strip())
        flush()
    return groups


def kmer_set(seq, k):
    seq = seq.upper()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else set()


def mash_identity(a, b, k):
    """Mash distance -> identity from k-mer Jaccard."""
    if not a or not b:
        return 0.0
    j = len(a & b) / len(a | b)
    if j <= 0:
        return 0.0
    if j >= 1:
        return 1.0
    d = -1.0 / k * math.log(2 * j / (1 + j))
    return max(0.0, 1.0 - d)


def group_stats(seqs, k):
    ksets = [kmer_set(s, k) for s, _ in seqs]
    n = len(ksets)
    if n <= 1:
        return n, 1.0, 1.0  # single copy: perfectly self-consistent
    ids = [mash_identity(a, b, k) for a, b in combinations(ksets, 2)]
    return n, sum(ids) / len(ids), min(ids)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--kmer", type=int, default=21)
    ap.add_argument("--threads", type=int, default=1)  # accepted for interface parity
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    groups = read_groups(args.fasta)
    stats = {}
    for gkey, seqs in groups.items():
        stats[gkey] = group_stats(seqs, args.kmer)

    with open(args.markers) as fh, open(args.out, "w") as out:
        header = fh.readline().rstrip("\n")
        out.write(header + "\tn_rep_copies\tmean_within_id\tmin_within_id\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rank, clade, cluster = f[0], f[1], f[2]
            n, mean_id, min_id = stats.get((rank, clade, cluster), (0, float("nan"), float("nan")))
            out.write(line.rstrip("\n") + f"\t{n}\t{mean_id:.4f}\t{min_id:.4f}\n")


if __name__ == "__main__":
    main()
