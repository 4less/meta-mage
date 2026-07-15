#!/usr/bin/env python3
"""
Build the genome manifest from the input samplesheet.

Input samplesheet (CSV or TSV, optional header) with three columns in any order:
    - accession        e.g. GCF_004353185.1_ASM435318v1_genomic
    - gtdb_lineage     d__..;p__..;c__..;o__..;f__..;g__..;s__..
    - path             path to the genome FASTA (.fna.gz)

Delimiter (tab vs comma) is auto-detected, and the three columns are identified
by content (lineage = the `;`/`__` field, path = the path-looking field), so the
column ORDER does not matter -- GTDB-downloader emits accession/path/lineage.

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

Optional --gtdb_metadata (GTDB bac120/ar53 metadata TSV, .gz ok) attaches each
genome's CheckM2 completeness/contamination. These feed the completeness-weighted
core-prevalence downstream: a genome missing a gene subtracts only its completeness
from the prevalence denominator, so a core gene absent merely because a MAG is
fragmentary is not demoted. Genomes with no metadata match default to completeness
1.0 (full weight = old behaviour).

Output TSV columns:
    idx  genome_id  is_rep  path  domain..species  completeness  contamination
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter

# Assembly-accession core, e.g. GCF_004353185.1 / GCA_000273195.1, extracted from
# whatever decoration the genome_id or GTDB accession carries.
_ACC_RE = re.compile(r"GC[AF]_\d+\.\d+")


def _acc_core(s):
    m = _ACC_RE.search(s)
    return m.group(0) if m else None


def load_gtdb_metadata(path):
    """accession-core -> (completeness_frac, contamination_frac).

    Keyed by both the RefSeq/GenBank forms (GCF_* from `accession`, GCA_* from
    `ncbi_genbank_assembly_accession`) so either style of genome_id resolves.
    """
    comp = {}
    if not path or os.path.basename(path) == "NO_FILE" or not os.path.exists(path):
        return comp
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        ci = col.get("checkm2_completeness", col.get("checkm_completeness"))
        ki = col.get("checkm2_contamination", col.get("checkm_contamination"))
        ai = col.get("accession")
        gi = col.get("ncbi_genbank_assembly_accession")
        if ci is None or ai is None:
            sys.stderr.write(f"GTDB metadata {path}: no completeness/accession "
                             "column found; skipping completeness join\n")
            return comp
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                c = float(f[ci]) / 100.0
                k = float(f[ki]) / 100.0 if ki is not None else 0.0
            except (ValueError, IndexError):
                continue
            for field in (f[ai], f[gi] if gi is not None and gi < len(f) else ""):
                core = _acc_core(field)
                if core:
                    comp[core] = (c, k)
    sys.stderr.write(f"GTDB metadata: loaded quality for {len(comp)} accessions\n")
    return comp

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
SPECIES_IDX = RANKS.index("species")


def parse_lineage(lineage_str):
    toks = [t.strip() for t in lineage_str.split(";")]
    return [toks[i] if i < len(toks) and toks[i] else "NA" for i in range(len(RANKS))]


def looks_like_header(row):
    joined = ",".join(row).lower()
    return "accession" in joined or "lineage" in joined or "path" in joined


def sniff_delimiter(path):
    """Tab if any tab appears in the first non-empty line, else comma."""
    with open(path, newline="") as fh:
        for line in fh:
            if line.strip():
                return "\t" if "\t" in line else ","
    return ","


PATH_SUFFIXES = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz", ".gz")


def classify_columns(fields):
    """Map fields to (accession, lineage, path, is_rep) regardless of order.

    lineage = the GTDB-lineage field (has ';' and a rank prefix like 'd__');
    path    = the field that looks like a filesystem path;
    is_rep  = a bare '0'/'1' flag column (GTDB-downloader's is_representative),
              or None when the sheet has no such column;
    accession = whatever remains. Returns None if a required column is missing.
    """
    li = next((i for i, f in enumerate(fields) if ";" in f and "__" in f), None)
    pi = next((i for i, f in enumerate(fields)
               if i != li and ("/" in f or f.lower().endswith(PATH_SUFFIXES))), None)
    ri = next((i for i, f in enumerate(fields)
               if i not in (li, pi) and f in ("0", "1")), None)
    ai = next((i for i in range(len(fields)) if i not in (li, pi, ri)), None)
    if None in (li, pi, ai):
        return None
    return fields[ai], fields[li], fields[pi], (fields[ri] if ri is not None else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome_sheet", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dropped",
                    help="TSV report of species dropped by the prefilter "
                         "(species, genus, n_genomes, min_required). Always "
                         "written (header only when nothing is dropped).")
    ap.add_argument("--min_genomes_per_species", type=int, default=10,
                    help="Drop genomes whose species has fewer than N genomes "
                         "in the input (prefilter). 1 = no filtering.")
    ap.add_argument("--gtdb_metadata",
                    help="GTDB metadata TSV(.gz) for CheckM2 completeness/"
                         "contamination. Omit (or NO_FILE) to weight every "
                         "genome as complete.")
    args = ap.parse_args()

    quality = load_gtdb_metadata(args.gtdb_metadata)

    delim = sniff_delimiter(args.genome_sheet)
    rows = []
    with open(args.genome_sheet, newline="") as fh:
        reader = csv.reader(fh, delimiter=delim)
        for i, row in enumerate(reader):
            row = [c.strip() for c in row if c.strip() != ""]
            if len(row) < 3:
                continue
            if i == 0 and looks_like_header(row):
                continue
            cols = classify_columns(row)
            if cols is None:
                sys.stderr.write(f"Skipping unparseable row {i}: {row}\n")
                continue
            rows.append(list(cols))

    if not rows:
        sys.exit(f"No genomes parsed from {args.genome_sheet}")

    # Prefilter: count genomes per species, then drop under-sampled species.
    ranks_by_row = [parse_lineage(lineage_str) for _, lineage_str, _, _ in rows]
    species_counts = Counter(r[SPECIES_IDX] for r in ranks_by_row)
    kept = [(row, ranks) for row, ranks in zip(rows, ranks_by_row)
            if species_counts[ranks[SPECIES_IDX]] >= args.min_genomes_per_species]

    n_dropped = len(rows) - len(kept)
    # species -> genus for every dropped (under-sampled) species, so the report
    # can break removals down per genus.
    GENUS_IDX = RANKS.index("genus")
    dropped_map = {
        r[SPECIES_IDX]: r[GENUS_IDX]
        for r in ranks_by_row
        if species_counts[r[SPECIES_IDX]] < args.min_genomes_per_species
    }
    dropped_species = sorted(dropped_map)
    if args.min_genomes_per_species > 1:
        sys.stderr.write(
            f"Prefilter --min_genomes_per_species={args.min_genomes_per_species}: "
            f"dropped {n_dropped} genome(s) across {len(dropped_species)} "
            f"under-sampled species: {', '.join(dropped_species) or '-'}\n"
        )
    # Always emit the drop report (even empty) so the downstream channel is stable.
    if args.dropped:
        with open(args.dropped, "w") as dr:
            dr.write("species\tgenus\tn_genomes\tmin_required\n")
            for sp in dropped_species:
                dr.write(f"{sp}\t{dropped_map[sp]}\t{species_counts[sp]}"
                         f"\t{args.min_genomes_per_species}\n")
    if not kept:
        sys.exit(
            f"No genomes left after --min_genomes_per_species="
            f"{args.min_genomes_per_species} prefilter"
        )

    n_matched = 0
    with open(args.out, "w") as out:
        out.write("\t".join(["idx", "genome_id", "is_rep", "path"] + RANKS
                            + ["completeness", "contamination"]) + "\n")
        for idx, ((accession, lineage_str, path, rep_flag), ranks) in enumerate(kept):
            # Prefer an explicit is_representative column; fall back to the
            # '.speciesrep' path/accession convention when the sheet lacks one.
            if rep_flag is not None:
                is_rep = rep_flag
            else:
                is_rep = "1" if (".speciesrep" in path or ".speciesrep" in accession) else "0"
            core = _acc_core(accession)
            comp, cont = quality.get(core, (1.0, 0.0)) if core else (1.0, 0.0)
            if quality and core in quality:
                n_matched += 1
            manifest_row = ([str(idx), accession, is_rep, os.path.abspath(path)]
                            + ranks + [f"{comp:.4f}", f"{cont:.4f}"])
            out.write("\t".join(manifest_row) + "\n")

    if quality:
        sys.stderr.write(f"Completeness join: matched {n_matched}/{len(kept)} "
                         "genomes to GTDB metadata (unmatched -> completeness 1.0)\n")
    sys.stderr.write(f"Wrote {len(kept)} genomes to {args.out}\n")


if __name__ == "__main__":
    main()
