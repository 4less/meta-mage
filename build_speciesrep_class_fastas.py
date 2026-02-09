#!/usr/bin/env python3
"""
Build pan-genome class FASTAs using only sequences from the speciesrep genome.

Outputs:
  - pan_genome_reference.core.speciesrep.fasta
  - pan_genome_reference.soft_core.speciesrep.fasta
  - pan_genome_reference.shell.speciesrep.fasta
  - pan_genome_reference.cloud.speciesrep.fasta

Usage:
    python3 build_speciesrep_class_fastas.py <panaroo_dir> <prokka_dir_or_speciesrep_ffn>
"""

import sys
from pathlib import Path
import pandas as pd


def read_presence_absence_table(panaroo_dir):
    candidates = [
        panaroo_dir / "gene_presence_absence.csv",
        panaroo_dir / "gene_presence_absence_roary.csv",
        panaroo_dir / "gene_presence_absence.tsv",
    ]
    for path in candidates:
        if path.exists():
            sep = "\t" if path.suffix == ".tsv" else ","
            return pd.read_csv(path, sep=sep, low_memory=False), path
    raise FileNotFoundError(
        "No gene_presence_absence file found. Checked: "
        "gene_presence_absence.csv, gene_presence_absence_roary.csv, gene_presence_absence.tsv"
    )


def find_speciesrep_column(df):
    rep_cols = [c for c in df.columns if ".speciesrep" in c]
    if not rep_cols:
        raise RuntimeError("No speciesrep column found in gene_presence_absence table")
    return rep_cols[0]


def find_speciesrep_ffn(prokka_dir):
    matches = sorted(prokka_dir.glob("*.speciesrep.ffn"))
    if not matches:
        raise FileNotFoundError(f"No '*.speciesrep.ffn' found in {prokka_dir}")
    if len(matches) > 1:
        raise RuntimeError(
            f"Multiple '*.speciesrep.ffn' found in {prokka_dir}: "
            + ", ".join(m.name for m in matches)
        )
    return matches[0]


def resolve_speciesrep_ffn(path_value):
    path = Path(path_value)
    if path.is_file():
        if not path.name.endswith(".ffn"):
            raise RuntimeError(f"Expected a .ffn file, got: {path}")
        return path
    if path.is_dir():
        return find_speciesrep_ffn(path)
    raise FileNotFoundError(f"{path} does not exist")


def read_rtab(rtab_file):
    if not rtab_file.exists():
        raise FileNotFoundError(f"{rtab_file} not found")
    return pd.read_csv(rtab_file, sep="\t", index_col=0)


def classify_genes(rtab):
    n = len(rtab.columns)
    if n == 0:
        return {k: set() for k in ("core", "soft_core", "shell", "cloud")}
    frac = (rtab == 1).sum(axis=1) / n
    return {
        "core": set(frac[frac == 1.0].index.tolist()),
        "soft_core": set(frac[(frac >= 0.95) & (frac < 1.0)].index.tolist()),
        "shell": set(frac[(frac >= 0.15) & (frac < 0.95)].index.tolist()),
        "cloud": set(frac[frac < 0.15].index.tolist()),
    }


def parse_fasta(fasta_path):
    seqs = {}
    cur_id = None
    cur_seq = []
    with open(fasta_path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            seqs[cur_id] = "".join(cur_seq)
    return seqs


def first_locus_tag(value):
    v = str(value).strip()
    if not v or v.lower() == "nan":
        return None
    return v.split(";")[0].strip()


def build_gene_to_locus(pa_df, rep_col):
    mapping = {}
    for _, row in pa_df.iterrows():
        gene = row["Gene"]
        tag = first_locus_tag(row[rep_col])
        if tag:
            mapping[gene] = tag
    return mapping


def write_class_fasta(out_file, genes, gene_to_locus, locus_to_seq):
    written = 0
    missing = 0
    with open(out_file, "w") as out:
        for gene in sorted(genes):
            locus = gene_to_locus.get(gene)
            if not locus:
                continue
            seq = locus_to_seq.get(locus)
            if not seq:
                missing += 1
                continue
            out.write(f">{gene}|{locus}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
            written += 1
    return written, missing


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    panaroo_dir = Path(sys.argv[1])
    rep_source = sys.argv[2]

    pa_df, pa_path = read_presence_absence_table(panaroo_dir)
    rep_col = find_speciesrep_column(pa_df)
    rtab = read_rtab(panaroo_dir / "gene_presence_absence.Rtab")
    classes = classify_genes(rtab)

    rep_ffn = resolve_speciesrep_ffn(rep_source)
    locus_to_seq = parse_fasta(rep_ffn)
    gene_to_locus = build_gene_to_locus(pa_df, rep_col)

    print(f"Presence/absence table: {pa_path}")
    print(f"Representative column: {rep_col}")
    print(f"Representative FFN: {rep_ffn}")

    for cls in ("core", "soft_core", "shell", "cloud"):
        out = panaroo_dir / f"pan_genome_reference.{cls}.speciesrep.fasta"
        written, missing = write_class_fasta(out, classes[cls], gene_to_locus, locus_to_seq)
        print(f"{cls}: wrote {written} sequences to {out}")
        if missing:
            print(f"{cls}: skipped {missing} genes (locus tag not found in speciesrep FFN)")


if __name__ == "__main__":
    main()
