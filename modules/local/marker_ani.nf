process MARKER_ANI {
    tag "marker_ani"
    label 'process_medium'
    publishDir "${params.outdir}/markers", mode: 'copy'
    cache 'lenient'

    // stdlib-python only (mash-style k-mer identity); runs in the base env like
    // DIVERSITY / REPORT, so no extra container is provisioned.

    input:
    path marker_fasta   // emitted marker CDS (nucleotide; rank|clade|cluster|genome)
    path scored         // markers table with the SCORE score column
    path spec_report    // specificity_report.tsv, or the NO_FILE placeholder

    output:
    path 'marker_ani.json', emit: ani

    script:
    // Pairwise marker-vs-marker ANI per species at three filtering stages, top
    // params.marker_ani_cap markers by score. Feeds the report's boxplot panels.
    """
    marker_ani.py \\
        --fasta ${marker_fasta} \\
        --scored ${scored} \\
        --specificity ${spec_report} \\
        --rank species \\
        --cap ${params.marker_ani_cap} \\
        --kmer ${params.marker_ani_kmer} \\
        --out marker_ani.json
    """

    stub:
    """
    printf '{"kmer":%s,"cap":%s,"stages":["specific-200","post-crossmap","post-recovery"],"species":{}}' ${params.marker_ani_kmer} ${params.marker_ani_cap} > marker_ani.json
    """
}
