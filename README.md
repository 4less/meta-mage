# meta-mage

Clade-specific marker gene discovery across GTDB, at scale (~800k genomes).

Given genome assemblies and their GTDB lineages, the pipeline finds genes that
are present in **> threshold** of an in-clade's genomes and in **< threshold** of
outside genomes, at **every taxonomic rank**, and emits a nucleotide
**classifier marker database**. It is a Nextflow (DSL2) pipeline built on nf-core
modules; the repository root *is* the pipeline.

## Input

A single samplesheet (`--genome_sheet`, CSV, header optional) with three columns:

```
accession,gtdb_lineage,path
GCF_004353185.1_...speciesrep,d__Bacteria;...;s__Clostridium cuniculi,/data/GCF_004353185.1_...fna.gz
```

A genome is treated as its species representative when `.speciesrep` appears in
its path or accession. Nucleotide marker sequences are drawn only from these reps.

## Running

```bash
nextflow run . --help
nextflow run . -profile test,docker                 # two-species fixture under data/
nextflow run . --genome_sheet sheet.csv -profile docker   # add ,slurm on a cluster
```

Wiring can be checked without the bio tools via `-stub-run`. The core selection
logic (`bin/*.py`) is stdlib-only and unit-testable directly.

## Design (why it scales)

| Stage | Module | Notes |
|-------|--------|-------|
| Manifest | `local/build_manifest` | Parses the samplesheet; assigns each genome a compact **integer id**; splits lineage into ranks. Single source of truth. |
| Gene calling | `nf-core/pyrodigal` | Proteins (`.faa`) for every genome + nucleotide CDS (`.fna`). The big **scatter**. |
| Reheader | `local/reheader` | Prefixes ids with `g<idx>_` so proteins/CDS are globally unique. |
| Vocabulary | `nf-core/mmseqs/{createdb,linclust,createtsv}` | One global **protein-space** clustering → gene-family vocabulary. Near-linear, disk-based. |
| Counts | `local/counts` | Presence/absence as **per-clade counts, never a matrix**. External `sort` groups by cluster; aggregates up the GTDB tree → all ranks in one pass. |
| Score | `local/score` | In/out prevalence per (cluster, rank, clade); threshold + cap per clade. |
| Payload | `local/emit_reps` | Nucleotide marker sequences pulled **only from in-clade species reps**. |
| Diversity | `local/diversity` | Within-clade nucleotide divergence per marker → #reps + threshold hints. |
| DB | `local/emit_db` | Final self-describing FASTA + marker table. |

**Key idea:** clustering stays in **protein space** (so homologs at higher ranks
stay one family), while the payload is **nucleotide, from species reps only**.
The only genome-scale structure held in memory is the compact `idx → lineage`
map; everything else streams or is externally sorted.

## Key parameters (see `nextflow.config`)

- `min_genomes_per_species` (10) — **prefilter**: drop genomes whose species has
  fewer than this many genomes before indexing (under-sampled species make
  within-species prevalence meaningless and inflate higher-rank clade sizes; 1 = off)
- `min_in_prevalence` (0.90), `max_out_prevalence` (0.02) — specificity thresholds
- `min_clade_size` (3), `max_markers_per_clade` (200)
- `linclust_min_seq_id` (0.5), `linclust_coverage` (0.8) — family granularity
- `div_kmer` (21) — within-clade distance k

`data/` holds the small two-species fixture used by `-profile test`. nf-core
modules are tracked in `modules.json`.
