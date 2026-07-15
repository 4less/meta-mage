process EMIT_REPS {
    tag "emit_reps"
    publishDir "${params.outdir}/markers", mode: 'copy'

    input:
    path markers           // rank clade cluster ...  (cluster ids tagged L:/S:)
    path clusters_loose    // rep <TAB> member  (loose clustering)
    path clusters_species  // rep <TAB> member  (species-tight clustering)
    path reps_ffn          // reheadered nucleotide CDS, species reps only
    path manifest

    output:
    path 'markers.nuc.fasta',  emit: marker_fasta
    path 'markers.emitted.tsv', emit: markers      // input markers minus too-short

    script:
    """
    emit_reps.py \\
        --markers ${markers} \\
        --clusters ${clusters_loose} \\
        --clusters_species ${clusters_species} \\
        --reps_ffn ${reps_ffn} \\
        --manifest ${manifest} \\
        --min_gene_len ${params.min_gene_len} \\
        --out_markers markers.emitted.tsv \\
        --out markers.nuc.fasta
    """
}
