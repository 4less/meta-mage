#!/usr/bin/env python3
"""
Emit nucleotide marker sequences, drawn ONLY from species-representative genomes.

For each selected marker (cluster) and its clade, collect the cluster members
that belong to species reps inside that clade, and write their CDS. That yields
one nucleotide sequence per constituent species rep -- a principled, GTDB-anchored
set of representatives per marker (captures within-clade diversity for free).

Cluster ids are namespaced by the clustering they came from ("L:" loose,
"S:" species; see aggregate_counts.py). A marker's members are looked up in the
matching clustering TSV, so species-rank markers resolve against the tight
clustering and higher ranks against the loose one.

Header: >{rank}|{clade}|{cluster}|{genome_id}
"""
import argparse

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def load_manifest(path):
    idx_is_rep = {}
    idx_gid = {}
    idx_lineage = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            idx = int(f[col["idx"]])
            idx_is_rep[idx] = f[col["is_rep"]] == "1"
            idx_gid[idx] = f[col["genome_id"]]
            idx_lineage[idx] = {r: f[col[r]] for r in RANKS}
    return idx_is_rep, idx_gid, idx_lineage


def load_markers(path):
    """tagged_cluster -> list of (rank, clade). Cluster id keeps its T:rep tag."""
    wanted = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rank, clade, cluster = f[0], f[1], f[2]
            wanted.setdefault(cluster, []).append((rank, clade))
    return wanted


def collect_members(clusters_path, tag, wanted, idx_is_rep, cluster_members):
    """Scan one rep<TAB>member TSV, recording rep-genome members of wanted
    clusters that carry this clustering's tag."""
    with open(clusters_path) as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")[:2]
            key = f"{tag}:{rep}"
            if key in wanted and key in cluster_members:
                g = genome_idx(member)
                if idx_is_rep.get(g, False):
                    cluster_members[key].append(member)


def genome_idx(protein_id):
    return int(protein_id[1:].split("_", 1)[0])


def load_fasta(path):
    seqs = {}
    cur = None
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers", required=True)
    ap.add_argument("--clusters", required=True,
                    help="loose clustering rep<TAB>member (tag L:)")
    ap.add_argument("--clusters_species", required=True,
                    help="species clustering rep<TAB>member (tag S:)")
    ap.add_argument("--reps_ffn", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    idx_is_rep, idx_gid, idx_lineage = load_manifest(args.manifest)
    wanted = load_markers(args.markers)
    # NOTE (scale): reps_ffn is reps-only and thus small vs the full corpus, but
    # at ~110k reps consider faidx random access instead of a full load.
    rep_seqs = load_fasta(args.reps_ffn)

    # Collect, for each wanted cluster, its rep member proteins. Each clustering
    # is scanned once and only contributes members to clusters carrying its tag.
    cluster_members = {c: [] for c in wanted}
    collect_members(args.clusters, "L", wanted, idx_is_rep, cluster_members)
    collect_members(args.clusters_species, "S", wanted, idx_is_rep, cluster_members)

    written = 0
    with open(args.out, "w") as out:
        for cluster, targets in wanted.items():
            members = cluster_members.get(cluster, [])
            for rank, clade in targets:
                for member in members:
                    g = genome_idx(member)
                    if idx_lineage[g][rank] != clade:
                        continue
                    seq = rep_seqs.get(member)
                    if not seq:
                        continue
                    out.write(f">{rank}|{clade}|{cluster}|{idx_gid[g]}\n")
                    for i in range(0, len(seq), 80):
                        out.write(seq[i:i + 80] + "\n")
                    written += 1
    print(f"emit_reps: wrote {written} marker sequences")


if __name__ == "__main__":
    main()
