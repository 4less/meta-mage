#!/usr/bin/env python3
"""
Turn protein cluster TSV(s) into per-clade presence counts, at every GTDB rank,
without ever materialising a genome x family matrix.

Rank-aware clustering: the vocabulary is defined at two protein identities.
  * The LOOSE clustering (e.g. 0.5) scaffolds family-and-above homology, where
    diverged orthologs must stay in one family.
  * The SPECIES clustering (e.g. 0.9) keeps paralogs / sister-species genes
    apart so species-rank specificity survives.
Each rank is counted from exactly one clustering: species from the tight pass,
every higher rank from the loose pass. Cluster ids are namespaced by their
source ("L:" loose, "S:" species) so downstream stages can route a marker back
to the clustering it came from.

Memory model:
  * idx -> lineage lives in RAM once (compact; clade strings interned).
  * clade sizes are plain counters (no per-clade genome sets); they depend only
    on the manifest, so they are identical across clusterings and written once.
  * Each cluster TSV is pre-sorted by cluster id, so we hold only the CURRENT
    cluster's set of genome indices while streaming.

For each cluster we emit, per rank, one row per clade that carries it:
    cluster  rank  clade  in_count  marker_total
where in_count = # in-clade genomes with the marker and marker_total = total
genomes with the marker. Out-of-clade prevalence is reconstructed downstream as
(marker_total - in_count) / (N - clade_size).
"""
import argparse
import sys
from sys import intern

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
SPECIES_RANK = "species"


def load_manifest(path):
    """idx -> (clade_at_each_rank tuple); plus clade sizes and N."""
    idx_lineage = {}
    clade_size = [dict() for _ in RANKS]  # one counter dict per rank
    n_total = 0
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            idx = int(f[col["idx"]])
            lineage = tuple(intern(f[col[r]]) for r in RANKS)
            idx_lineage[idx] = lineage
            n_total += 1
            for ri, clade in enumerate(lineage):
                clade_size[ri][clade] = clade_size[ri].get(clade, 0) + 1
    return idx_lineage, clade_size, n_total


def genome_idx_from_protein(protein_id):
    # g<idx>_<original token>
    return int(protein_id[1:].split("_", 1)[0])


def flush_cluster(cluster, genomes, idx_lineage, rank_indices, tag, out):
    """genomes: set of genome idx carrying this cluster.

    rank_indices: which RANKS positions to emit rows for (this clustering owns
    those ranks). tag: cluster-id namespace prefix ("L"/"S")."""
    marker_total = len(genomes)
    if marker_total == 0:
        return
    for ri in rank_indices:
        rank = RANKS[ri]
        clade_counts = {}
        for g in genomes:
            clade = idx_lineage[g][ri]
            clade_counts[clade] = clade_counts.get(clade, 0) + 1
        for clade, in_count in clade_counts.items():
            if clade == "NA":
                continue
            out.write(f"{tag}:{cluster}\t{rank}\t{clade}\t{in_count}\t{marker_total}\n")


def stream_clustering(clusters_path, idx_lineage, rank_indices, tag, out):
    """Stream one sorted rep<TAB>member TSV, emitting counts for rank_indices."""
    if not rank_indices:
        return
    with open(clusters_path) as fh:
        cur = None
        genomes = set()
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")[:2]
            if rep != cur:
                if cur is not None:
                    flush_cluster(cur, genomes, idx_lineage, rank_indices, tag, out)
                cur = rep
                genomes = set()
            try:
                genomes.add(genome_idx_from_protein(member))
            except (ValueError, IndexError):
                sys.stderr.write(f"skip unparseable member: {member}\n")
        if cur is not None:
            flush_cluster(cur, genomes, idx_lineage, rank_indices, tag, out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True,
                    help="loose clustering rep<TAB>member, sorted by rep")
    ap.add_argument("--clusters_species",
                    help="species clustering rep<TAB>member, sorted by rep. "
                         "If omitted, the loose clustering owns every rank.")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--clade_sizes", required=True)
    args = ap.parse_args()

    idx_lineage, clade_size, n_total = load_manifest(args.manifest)

    # Clade sizes come from the manifest alone -> identical across clusterings.
    with open(args.clade_sizes, "w") as cs:
        cs.write("rank\tclade\tsize\n")
        cs.write(f"__TOTAL__\t__ALL__\t{n_total}\n")
        for ri, rank in enumerate(RANKS):
            for clade, size in clade_size[ri].items():
                if clade == "NA":
                    continue
                cs.write(f"{rank}\t{clade}\t{size}\n")

    species_ri = RANKS.index(SPECIES_RANK)
    if args.clusters_species:
        loose_ranks = [ri for ri in range(len(RANKS)) if ri != species_ri]
        species_ranks = [species_ri]
    else:
        loose_ranks = list(range(len(RANKS)))
        species_ranks = []

    with open(args.counts, "w") as out:
        out.write("cluster\trank\tclade\tin_count\tmarker_total\n")
        stream_clustering(args.clusters, idx_lineage, loose_ranks, "L", out)
        if species_ranks:
            stream_clustering(args.clusters_species, idx_lineage,
                              species_ranks, "S", out)


if __name__ == "__main__":
    main()
