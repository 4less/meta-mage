#!/usr/bin/env python3
"""
Emit nucleotide marker sequences, one per constituent species of the marker's clade.

For each selected marker (cluster) and its clade, collect the cluster members and,
for each species in the clade, emit the CDS from the BEST genome carrying the gene:
the species representative if it has the marker, otherwise the highest-CheckM2-
completeness genome that does (length ties broken toward the rep, then completeness,
then idx for determinism). This decouples marker emission from a single rep's gene
content -- a near-universal gene absent only because the rep genome is fragmentary
is no longer silently dropped (the rep-completeness failure mode).

Sequences are drawn from all_cds (every genome's CDS, reheadered g<idx>_<orig>).
Only the member proteins of selected clusters are held in memory: a single filtered
streaming pass over all_cds keeps just those sequences, so peak RAM tracks the
marker payload, not the whole corpus.

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
    idx_comp = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        ci = col.get("completeness")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            idx = int(f[col["idx"]])
            idx_is_rep[idx] = f[col["is_rep"]] == "1"
            idx_gid[idx] = f[col["genome_id"]]
            idx_lineage[idx] = {r: f[col[r]] for r in RANKS}
            idx_comp[idx] = (float(f[ci]) if ci is not None and f[ci] != "" else 1.0)
    return idx_is_rep, idx_gid, idx_lineage, idx_comp


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


def collect_members(clusters_path, tag, cluster_members):
    """Scan one rep<TAB>member TSV, recording ALL members of wanted clusters that
    carry this clustering's tag (not just species reps -- the best-genome fallback
    needs every genome that carries the gene as a candidate)."""
    with open(clusters_path) as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")[:2]
            key = f"{tag}:{rep}"
            if key in cluster_members:
                cluster_members[key].append(member)


def genome_idx(protein_id):
    return int(protein_id[1:].split("_", 1)[0])


def load_seqs_subset(path, needed):
    """Single streaming pass over a (large) CDS FASTA, keeping only sequences whose
    id is in `needed`. Memory is bounded by the marker payload, not the whole corpus.
    NOTE (scale): a faidx random-access read of just `needed` avoids even this one
    linear scan of all_cds -- worth it once all_cds is 100s of GB."""
    seqs = {}
    cur = None
    keep = False
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None and keep:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                keep = cur in needed
                buf = []
            elif keep:
                buf.append(line.strip())
        if cur is not None and keep:
            seqs[cur] = "".join(buf)
    return seqs


def genome_rank_key(g, idx_is_rep, idx_comp):
    """Best genome first: rep, then highest completeness, then lowest idx."""
    return (0 if idx_is_rep.get(g, False) else 1, -idx_comp.get(g, 0.0), g)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers", required=True)
    ap.add_argument("--clusters", required=True,
                    help="loose clustering rep<TAB>member (tag L:)")
    ap.add_argument("--clusters_species", required=True,
                    help="species clustering rep<TAB>member (tag S:)")
    ap.add_argument("--all_cds", required=True,
                    help="reheadered nucleotide CDS for EVERY genome (g<idx>_<orig>)")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out_markers",
                    help="filtered marker table: input markers minus those left "
                         "with no sequence >= min_gene_len (keeps the DB consistent)")
    ap.add_argument("--min_gene_len", type=int, default=0,
                    help="drop marker CDS shorter than this many bp (0 = off)")
    args = ap.parse_args()

    idx_is_rep, idx_gid, idx_lineage, idx_comp = load_manifest(args.manifest)
    wanted = load_markers(args.markers)

    # Collect, for each wanted cluster, ALL its member proteins. Each clustering is
    # scanned once and only contributes members to clusters carrying its tag.
    cluster_members = {c: [] for c in wanted}
    collect_members(args.clusters, "L", cluster_members)
    collect_members(args.clusters_species, "S", cluster_members)

    # Pull only the member sequences of selected clusters out of all_cds.
    needed = set()
    for members in cluster_members.values():
        needed.update(members)
    seqs = load_seqs_subset(args.all_cds, needed)

    written = 0
    short = 0
    fallback = 0              # species whose CDS came from a non-rep genome
    emitted_markers = set()   # (rank, clade, cluster) with >= 1 kept sequence
    with open(args.out, "w") as out:
        for cluster, targets in wanted.items():
            members = cluster_members.get(cluster, [])
            for rank, clade in targets:
                # group member proteins by species -> genome, within this clade
                by_species = {}
                for m in members:
                    g = genome_idx(m)
                    if idx_lineage[g][rank] != clade:
                        continue
                    sp = idx_lineage[g]["species"]
                    by_species.setdefault(sp, {}).setdefault(g, []).append(m)
                # one sequence per constituent species, from the best genome that
                # actually carries a long-enough CDS (rep preferred, else fall back)
                for genomes in by_species.values():
                    for g in sorted(genomes,
                                    key=lambda x: genome_rank_key(x, idx_is_rep, idx_comp)):
                        emitted_any = False
                        for m in genomes[g]:
                            seq = seqs.get(m)
                            if not seq:
                                continue
                            if args.min_gene_len and len(seq) < args.min_gene_len:
                                short += 1
                                continue
                            out.write(f">{rank}|{clade}|{cluster}|{idx_gid[g]}\n")
                            for i in range(0, len(seq), 80):
                                out.write(seq[i:i + 80] + "\n")
                            written += 1
                            emitted_any = True
                            emitted_markers.add((rank, clade, cluster))
                        if emitted_any:
                            if not idx_is_rep.get(g, False):
                                fallback += 1
                            break   # this species is covered; skip lesser genomes

    # Filtered marker table: keep only markers that retained a sequence. A marker
    # whose every candidate CDS was below min_gene_len is dropped from the DB.
    dropped_markers = 0
    if args.out_markers:
        with open(args.markers) as fh, open(args.out_markers, "w") as mout:
            mout.write(fh.readline())  # header
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if (f[0], f[1], f[2]) in emitted_markers:
                    mout.write(line)
                else:
                    dropped_markers += 1
    print(f"emit_reps: wrote {written} marker sequences "
          f"({short} dropped < {args.min_gene_len}bp; "
          f"{fallback} species drawn from a non-rep genome; "
          f"{dropped_markers} markers dropped for having no long-enough sequence)")


if __name__ == "__main__":
    main()
