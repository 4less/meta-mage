#!/usr/bin/env python3
"""
Build the genome manifest from the input samplesheet.

Input samplesheet (CSV, optional header) with three columns:
    1. accession       e.g. GCF_004353185.1_ASM435318v1_genomic
    2. gtdb_lineage    d__..;p__..;c__..;o__..;f__..;g__..;s__..
    3. path            path to the genome FASTA (.fna.gz)

Assigns each genome a compact integer index (everything downstream keys on it,
memory-efficient at 800k scale) and splits the lineage into ranks.

A genome is treated as its species' representative if '.speciesrep' appears in
its path or accession (the GTDB-download convention used here).

Prefilter: `--min_genomes_per_species N` drops every genome whose species is
represented by fewer than N genomes in the input, *before* indexing. Under-sampled
species make within-species prevalence meaningless (with 1 genome it is trivially
0 or 1), inflate higher-rank clade sizes, and pollute the out-clade denominator.
Filtering here means the integer ids stay compact and no downstream stage ever
sees the dropped genomes. The species-representative of a dropped species goes
too. Default 1 = no filtering.

Output TSV columns:
    idx  genome_id  is_rep  path  domain phylum class order family genus species
"""
import argparse
import csv
import os
import sys
from collections import Counter

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
SPECIES_IDX = RANKS.index("species")


def parse_lineage(lineage_str):
    toks = [t.strip() for t in lineage_str.split(";")]
    return [toks[i] if i < len(toks) and toks[i] else "NA" for i in range(len(RANKS))]


def looks_like_header(row):
    joined = ",".join(row).lower()
    return "accession" in joined or "lineage" in joined or "path" in joined


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome_sheet", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min_genomes_per_species", type=int, default=10,
                    help="Drop genomes whose species has fewer than N genomes "
                         "in the input (prefilter). 1 = no filtering.")
    args = ap.parse_args()

    rows = []
    with open(args.genome_sheet, newline="") as fh:
        reader = csv.reader(fh)
        for i, row in enumerate(reader):
            row = [c.strip() for c in row if c.strip() != ""]
            if len(row) < 3:
                continue
            if i == 0 and looks_like_header(row):
                continue
            rows.append(row[:3])

    if not rows:
        sys.exit(f"No genomes parsed from {args.genome_sheet}")

    # Prefilter: count genomes per species, then drop under-sampled species.
    ranks_by_row = [parse_lineage(lineage_str) for _, lineage_str, _ in rows]
    species_counts = Counter(r[SPECIES_IDX] for r in ranks_by_row)
    kept = [(row, ranks) for row, ranks in zip(rows, ranks_by_row)
            if species_counts[ranks[SPECIES_IDX]] >= args.min_genomes_per_species]

    n_dropped = len(rows) - len(kept)
    if args.min_genomes_per_species > 1:
        dropped_species = sorted(
            {r[SPECIES_IDX] for _, r in zip(rows, ranks_by_row)
             if species_counts[r[SPECIES_IDX]] < args.min_genomes_per_species}
        )
        sys.stderr.write(
            f"Prefilter --min_genomes_per_species={args.min_genomes_per_species}: "
            f"dropped {n_dropped} genome(s) across {len(dropped_species)} "
            f"under-sampled species: {', '.join(dropped_species) or '-'}\n"
        )
    if not kept:
        sys.exit(
            f"No genomes left after --min_genomes_per_species="
            f"{args.min_genomes_per_species} prefilter"
        )

    with open(args.out, "w") as out:
        out.write("\t".join(["idx", "genome_id", "is_rep", "path"] + RANKS) + "\n")
        for idx, ((accession, lineage_str, path), ranks) in enumerate(kept):
            is_rep = "1" if (".speciesrep" in path or ".speciesrep" in accession) else "0"
            manifest_row = [str(idx), accession, is_rep, os.path.abspath(path)] + ranks
            out.write("\t".join(manifest_row) + "\n")

    sys.stderr.write(f"Wrote {len(kept)} genomes to {args.out}\n")


if __name__ == "__main__":
    main()
