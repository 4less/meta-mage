# Core Gene Extraction Script

## Overview

`extract_core_genes.py` processes panaroo output to identify and extract core genes—genes present in all genomes and in the representative genome.

## How It Works

The script implements a two-stage filtering approach:

1. **Core gene identification**: Genes labeled as "core" with presence value `1` in ALL genomes (from `gene_presence_absence.Rtab`)
2. **Representative genome validation**: Further filters to ensure the gene is present in the representative genome (identified by `.speciesrep.fna.gz` suffix)

## Input Files

The script expects panaroo output in the following structure:
```
panaroo_output/
├── gene_presence_absence.tsv    # Gene metadata and annotations
└── gene_presence_absence.Rtab   # Binary presence/absence matrix (0/1)
```

### gene_presence_absence.Rtab Format
```
Gene                    genome1.fna.gz    genome2.fna.gz    rep.speciesrep.fna.gz    ...
group_1758              1                 1                 1                        ...
group_1756              1                 1                 1                        ...
```

## Usage

### Basic Usage
```bash
python3 extract_core_genes.py <panaroo_dir> <output_file>
```

### Example
```bash
python3 extract_core_genes.py ./output/s__Clostridium_cuniculi/panaroo ./core_genes.tsv
```

### From Workflow
The script is automatically run after panaroo in `workflow.sh`:
```bash
python3 extract_core_genes.py "$PANA" "$BASE_OUT/core_genes.tsv"
```

## Output

Generates a TSV file with core genes and their metadata:

```
Gene            Non-unique Gene name    Annotation
group_1758                              hypothetical protein
gerBA           gerBA                   Spore germination protein B1
pncC            pncC                    Nicotinamide-nucleotide amidohydrolase PncC
fusA            fusA                    Elongation factor G
...
```

## Algorithm

1. Read `gene_presence_absence.Rtab` to identify the representative genome
2. Filter for rows where ALL columns (genomes) have value `1` (core genes)
3. Verify the representative genome has the gene
4. Extract gene annotations from `gene_presence_absence.tsv`
5. Output deduplicated, sorted list of core genes

## Dependencies

- pandas
- Python 3.6+

Install pandas if needed:
```bash
pip install pandas
```

## Notes

- The representative genome is identified by the `.speciesrep.fna.gz` suffix in column names
- If no representative genome is found, a warning is printed but processing continues
- Core genes are deduplicated and sorted alphabetically
