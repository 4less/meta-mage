process SPECIFICITY {
    tag "specificity"
    publishDir "${params.outdir}/markers", mode: 'copy'

    input:
    path hits          // crossmap.m8: query target pident alnlen
    path idmap         // query.map: m<N> <TAB> original header
    path markers       // scored marker table
    path marker_fasta  // emitted marker CDS
    path manifest

    output:
    path 'markers.specific.tsv',   emit: markers
    path 'markers.specific.fasta', emit: marker_fasta
    path 'specificity_report.tsv', emit: report

    script:
    """
    specificity_guard.py \\
        --hits ${hits} \\
        --idmap ${idmap} \\
        --markers ${markers} \\
        --fasta ${marker_fasta} \\
        --manifest ${manifest} \\
        --min_id ${params.specificity_min_id} \\
        --min_aln ${params.specificity_min_aln} \\
        --out_markers markers.specific.tsv \\
        --out_fasta markers.specific.fasta \\
        --report specificity_report.tsv
    """
}
