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
    cluster  rank  clade  in_count  marker_total  in_score
where in_count = # in-clade genomes with the marker, marker_total = total genomes
with the marker, and in_score = sum of those in-clade genomes' completeness scores.
Out-of-clade prevalence is reconstructed downstream as
(marker_total - in_count) / (N - clade_size).

clade_sizes.tsv carries, per clade, both the genome count (size) and score_sum =
sum of completeness over the clade's genomes. The completeness-weighted core
prevalence downstream is (size - (score_sum - in_score)) / size: present genomes
count fully, each MISSING genome subtracts only its completeness. With no metadata
(every completeness = 1.0) score_sum == size and in_score == in_count, so the
weighted value collapses to the plain in_count / size.
"""
import argparse
import sys
from sys import intern

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
SPECIES_RANK = "species"


def load_manifest(path):
    """idx -> lineage tuple, idx -> completeness score; plus clade sizes/score_sums."""
    idx_lineage = {}
    idx_score = {}
    clade_size = [dict() for _ in RANKS]        # one counter dict per rank
    clade_score = [dict() for _ in RANKS]       # parallel sum of completeness
    n_total = 0
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        has_comp = "completeness" in col
        for line in fh:
            f = line.rstrip("\n").split("\t")
            idx = int(f[col["idx"]])
            lineage = tuple(intern(f[col[r]]) for r in RANKS)
            idx_lineage[idx] = lineage
            score = float(f[col["completeness"]]) if has_comp else 1.0
            idx_score[idx] = score
            n_total += 1
            for ri, clade in enumerate(lineage):
                clade_size[ri][clade] = clade_size[ri].get(clade, 0) + 1
                clade_score[ri][clade] = clade_score[ri].get(clade, 0.0) + score
    return idx_lineage, idx_score, clade_size, clade_score, n_total


def genome_idx_from_protein(protein_id):
    # g<idx>_<original token>
    return int(protein_id[1:].split("_", 1)[0])


def flush_cluster(cluster, genomes, idx_lineage, idx_score, rank_indices, tag, out):
    """genomes: set of genome idx carrying this cluster.

    rank_indices: which RANKS positions to emit rows for (this clustering owns
    those ranks). tag: cluster-id namespace prefix ("L"/"S"). Emits in_score =
    sum of present in-clade genomes' completeness alongside in_count."""
    marker_total = len(genomes)
    if marker_total == 0:
        return
    for ri in rank_indices:
        rank = RANKS[ri]
        clade_counts = {}
        clade_scores = {}
        for g in genomes:
            clade = idx_lineage[g][ri]
            clade_counts[clade] = clade_counts.get(clade, 0) + 1
            clade_scores[clade] = clade_scores.get(clade, 0.0) + idx_score[g]
        for clade, in_count in clade_counts.items():
            if clade == "NA":
                continue
            out.write(f"{tag}:{cluster}\t{rank}\t{clade}\t{in_count}\t"
                      f"{marker_total}\t{clade_scores[clade]:.4f}\n")


def stream_clustering(clusters_path, idx_lineage, idx_score, rank_indices, tag, out):
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
                    flush_cluster(cur, genomes, idx_lineage, idx_score,
                                  rank_indices, tag, out)
                cur = rep
                genomes = set()
            try:
                genomes.add(genome_idx_from_protein(member))
            except (ValueError, IndexError):
                sys.stderr.write(f"skip unparseable member: {member}\n")
        if cur is not None:
            flush_cluster(cur, genomes, idx_lineage, idx_score,
                          rank_indices, tag, out)


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

    idx_lineage, idx_score, clade_size, clade_score, n_total = \
        load_manifest(args.manifest)
    total_score = sum(clade_score[RANKS.index("domain")].values())

    # Clade sizes come from the manifest alone -> identical across clusterings.
    # score_sum = sum of completeness over the clade's genomes (== size when no
    # metadata was joined).
    with open(args.clade_sizes, "w") as cs:
        cs.write("rank\tclade\tsize\tscore_sum\n")
        cs.write(f"__TOTAL__\t__ALL__\t{n_total}\t{total_score:.4f}\n")
        for ri, rank in enumerate(RANKS):
            for clade, size in clade_size[ri].items():
                if clade == "NA":
                    continue
                cs.write(f"{rank}\t{clade}\t{size}\t{clade_score[ri][clade]:.4f}\n")

    species_ri = RANKS.index(SPECIES_RANK)
    if args.clusters_species:
        loose_ranks = [ri for ri in range(len(RANKS)) if ri != species_ri]
        species_ranks = [species_ri]
    else:
        loose_ranks = list(range(len(RANKS)))
        species_ranks = []

    with open(args.counts, "w") as out:
        out.write("cluster\trank\tclade\tin_count\tmarker_total\tin_score\n")
        stream_clustering(args.clusters, idx_lineage, idx_score,
                          loose_ranks, "L", out)
        if species_ranks:
            stream_clustering(args.clusters_species, idx_lineage, idx_score,
                              species_ranks, "S", out)


if __name__ == "__main__":
    main()
