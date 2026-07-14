#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MARKERS } from './workflows/markers.nf'

def helpMessage() {
    log.info """
    ============================================================
     meta-mage/markers  —  clade-specific marker gene discovery
    ============================================================

    Usage:
      nextflow run . --genome_sheet <csv> [options] -profile docker

    Required:
      --genome_sheet  CSV with three columns (header optional):
                        accession, gtdb_lineage (d__..;..;s__..), path-to-genome
                      A genome is the species rep if '.speciesrep' is in its path.

    Options (defaults in nextflow.config):
      --outdir                 Results directory                 [${params.outdir}]
      --min_in_prevalence      Min in-clade prevalence           [${params.min_in_prevalence}]
      --max_out_prevalence     Max outside-clade prevalence      [${params.max_out_prevalence}]
      --min_clade_size         Skip clades smaller than this     [${params.min_clade_size}]
      --max_markers_per_clade  Cap markers per clade             [${params.max_markers_per_clade}]
      --core_prevalence        "Core gene" cutoff for the report [${params.core_prevalence}]
      --linclust_min_seq_id    Protein cluster identity (loose)  [${params.linclust_min_seq_id}]
      --linclust_min_seq_id_species  Species-rank cluster id     [${params.linclust_min_seq_id_species}]
      --linclust_coverage      Protein cluster coverage          [${params.linclust_coverage}]
      --div_kmer               k for within-clade distance       [${params.div_kmer}]
      --specificity            Nucleotide cross-map guard on/off  [${params.specificity}]
      --specificity_min_id     Off-target hit identity            [${params.specificity_min_id}]
      --specificity_min_aln    Off-target hit min aln length (bp) [${params.specificity_min_aln}]

    Profiles:
      -profile test            Two-species fixture under data/
      -profile docker | conda  Provision the bio tools (pyrodigal, mmseqs)
      -profile slurm           Submit to a SLURM cluster

    Other:
      -stub-run                Validate wiring without running the tools
      --help                   Show this message
    """.stripIndent()
}

workflow {
    if( params.help ) {
        helpMessage()
        return
    }

    if( !params.genome_sheet )
        error "Provide --genome_sheet (CSV: accession,gtdb_lineage,path), or use -profile test. See --help."

    MARKERS()
}
