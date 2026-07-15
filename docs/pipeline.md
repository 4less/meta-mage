# meta-mage — pipeline flow

How the `markers` workflow ([workflows/markers.nf](../workflows/markers.nf)) turns genome
assemblies + GTDB lineages into a nucleotide **classifier marker database**.

Legend: **rounded** = process/tool · **hexagon** = filter/selector (drops or
subsets data) · **parallelogram** = data artifact · **green** = published output ·
dashed edges = the optional nucleotide specificity guard (`--specificity`).

A PDF of the *actual* run DAG is emitted every run at
`${outdir}/pipeline_info/pipeline_dag.pdf` (Nextflow + graphviz); the diagram below
is the curated illustration.

```mermaid
flowchart TD
    %% ---------- Inputs ----------
    IN[/"genome_sheet CSV<br/>accession, gtdb_lineage, path"/]:::data
    GTM[/"GTDB metadata TSV (optional)<br/>checkm2_completeness"/]:::data

    %% ---------- Manifest ----------
    IN --> BM("BUILD_MANIFEST<br/>build_manifest.py")
    GTM --> BM
    PREF{{"prefilter: drop species with<br/>&lt; min_genomes_per_species (10)<br/>assign idx · split ranks · join completeness"}}:::filt
    BM --> PREF --> MAN[/"manifest.tsv<br/>idx · genome_id · is_rep · lineage · completeness"/]:::data

    %% ---------- Per-genome scatter ----------
    MAN -->|"splitCsv → per-genome (SCATTER)"| PYR("PYRODIGAL  -p meta<br/>gene calling")
    PYR --> FAA[/".faa proteins"/]:::data
    PYR --> FNA[/".fna nucleotide CDS"/]:::data

    %% ---------- Protein vocabulary ----------
    FAA --> RHF("REHEADER_FAA<br/>prefix g&lt;idx&gt;_ → globally unique")
    RHF -->|collectFile| ALLFAA[/"all_proteins.faa"/]:::data
    ALLFAA --> DB("MMSEQS_CREATEDB<br/>protein DB")

    DB --> LL("MMSEQS_LINCLUST_LOOSE<br/>--min-seq-id 0.5 -c 0.8")
    DB --> LS("MMSEQS_LINCLUST_SPECIES<br/>tight id, UniRef90-like")
    LL --> CTL("CREATETSV_LOOSE") --> CLU[/"clusters_loose<br/>rep⇥member · owns ranks &gt; species"/]:::data
    LS --> CTS("CREATETSV_SPECIES") --> CLS[/"clusters_species<br/>rep⇥member · owns species rank"/]:::data

    %% ---------- Counts + scoring ----------
    CLU --> CNT("COUNTS<br/>LC_ALL=C sort + aggregate_counts.py")
    CLS --> CNT
    MAN --> CNT
    CNT --> CO[/"counts.tsv<br/>cluster·rank·clade·in_count·in_score"/]:::data
    CNT --> CS[/"clade_sizes.tsv<br/>size · score_sum (Σ completeness)"/]:::data

    CO --> SC("SCORE<br/>score_markers.py")
    CS --> SC
    SEL{{"select per (cluster,rank,clade):<br/>completeness-weighted in-prev ≥ 0.80<br/>out-prev ≤ 0.02 · clade_size ≥ 3<br/>rank by in·(1-out) · cap 100 / clade"}}:::filt
    SC --> SEL --> MK[/"markers.tsv"/]:::data

    %% ---------- Nucleotide payload ----------
    FNA --> RHN("REHEADER_FNA<br/>prefix g&lt;idx&gt;_")
    RHN --> REPFILT{{"filter is_rep:<br/>species reps only"}}:::filt
    REPFILT -->|collectFile| REPS[/"reps.ffn"/]:::data
    RHN -->|collectFile| ALLCDS[/"all_cds.ffn<br/>every genome"/]:::data

    MK --> ER("EMIT_REPS<br/>emit_reps.py")
    CLU --> ER
    CLS --> ER
    REPS --> ER
    MAN --> ER
    LENF{{"drop marker CDS<br/>&lt; min_gene_len (300 bp)"}}:::filt
    ER --> LENF
    LENF --> MFA[/"markers.nuc.fasta"/]:::data
    LENF --> MKE[/"markers.emitted.tsv<br/>length-filtered marker set"/]:::data

    %% ---------- Optional specificity guard ----------
    MFA -.-> CM("CROSSMAP<br/>mmseqs easy-search · nuc-vs-nuc (type 3)")
    ALLCDS -.-> CM
    CM -.-> HITS[/"crossmap.m8 + query.map"/]:::data
    HITS -.-> SP("SPECIFICITY<br/>specificity_guard.py")
    MKE -.-> SP
    MFA -.-> SP
    MAN -.-> SP
    SPFILT{{"drop markers cross-mapping off-target:<br/>pident ≥ specificity_min_id<br/>alnlen ≥ specificity_min_aln"}}:::filt
    SP -.-> SPFILT
    SPFILT -.-> SPOUT[/"markers.specific.tsv/.fasta"/]:::data

    %% ---------- Diversity + DB ----------
    MKE ==>|"--specificity off"| DIV
    MFA ==>|"--specificity off"| DIV
    SPOUT --> DIV("DIVERSITY<br/>diversity.py · mash k=21")
    DIV --> MDT[/"markers.diversity.tsv"/]:::data
    MDT --> ED("EMIT_DB<br/>emit_db.py · dedup + annotate")
    MFA --> ED
    ED --> DBOUT[/"marker_db.fasta + marker_db.tsv"/]:::out

    %% ---------- Report ----------
    MAN --> RPT("REPORT<br/>build_report.py")
    CO --> RPT
    CS --> RPT
    MKE --> RPT
    SPOUT -.-> RPT
    RPT --> RPTO[/"report.html"/]:::out

    %% ---------- Low-marker diagnostics ----------
    subgraph LMD ["low-marker diagnostics · species with &lt; low_marker_threshold (50) markers"]
      direction TB
      LMM("LOW_MARKER_MASKING<br/>low_marker_masking.py")
      MSKR("MASKING_REPORT<br/>masking_report.py")
      ASK("ANI_SKANI<br/>skani triangle (container)")
      AG("ANI_GAP<br/>ani_gap.py")
      LMR("LOW_MARKER_REPORT<br/>low_marker_report.py")
      LMM --> MSKR
      LMM --> AG
      ASK --> AG
      AG --> LMR
      LMM --> LMR
    end
    SPOUT -.-> LMM
    HITS -.-> LMM
    REPS -.-> LMM
    MAN -.-> LMM
    MSKR --> MSKRO[/"masking_report.html"/]:::out
    LMR --> LMRO[/"low_marker_report.html"/]:::out

    classDef data fill:#eef4ff,stroke:#5b7fb5,color:#12325e;
    classDef filt fill:#fff2e0,stroke:#d68a2b,color:#5e3a0f;
    classDef out  fill:#e6f6ea,stroke:#3f9d5a,color:#12401f;
```

## Key idea

Clustering stays in **protein space** (homologs at higher ranks stay one family),
while the marker payload is **nucleotide, from species reps only**. Two clusterings
run over the same protein DB: a **loose** one scaffolds family-and-above, a **tight**
one owns the species rank so sister-species/paralogs stay apart. The only
genome-scale structure held in memory is the compact `idx → lineage` map; everything
else streams or is externally sorted (the 800k-genome-safe path).

**Completeness weighting** (when `--gtdb_metadata` is set): core prevalence is
`(N − Σ completeness of the genomes missing the gene) / N`, so a core gene absent
only because a genome is fragmentary is not demoted. With no metadata every genome
weighs 1.0 and this reduces to plain presence/N.

**Low-marker diagnostics**: any species that ends below `low_marker_threshold`
markers is examined for (1) whether its cross-map-dropped markers could be salvaged
by masking, and (2) whether it has a within-vs-between ANI gap in its genus (a
re-merge signal).

See [README.md](../README.md) for the stage table and parameter reference.
