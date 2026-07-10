process EMIT_DB {
    tag "emit_db"
    publishDir "${params.outdir}/db", mode: 'copy'

    input:
    path marker_table   // markers.diversity.tsv (scored + diversity columns)
    path marker_fasta   // markers.nuc.fasta

    output:
    path 'marker_db.fasta',  emit: db
    path 'marker_db.tsv',    emit: table

    script:
    // Finalise the classifier DB: annotate headers with prevalence + diversity,
    // drop exact-duplicate copies within a marker.
    """
    emit_db.py \\
        --table ${marker_table} \\
        --fasta ${marker_fasta} \\
        --out_fasta marker_db.fasta \\
        --out_table marker_db.tsv
    """
}
