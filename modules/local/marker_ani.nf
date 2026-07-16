process MARKER_ANI {
    tag "marker_ani"
    publishDir "${params.outdir}/markers", mode: 'copy'
    cache 'lenient'

    // stdlib-python only (mash-style k-mer distance); runs in the base env like
    // DIVERSITY / REPORT, so no container is provisioned. Streams the big all-CDS
    // and species-cluster files once (see conf/base.config for its resources).

    input:
    path scored             // markers table with the SCORE score column
    path spec_report        // specificity_report.tsv, or the NO_FILE placeholder
    path clusters_species   // clusters_species.tsv: rep <TAB> member
    path manifest           // manifest.tsv (idx -> species)
    path all_cds            // all-genome CDS FASTA (nucleotide, >g<idx>_ headers)

    output:
    path 'marker_ani.json', emit: ani

    script:
    // Per-species WITHIN-marker nucleotide distance: for each top-scoring marker,
    // pairwise distance among its copies across the species' genomes, at three
    // filtering stages. Feeds the report's boxplot panels.
    """
    marker_ani.py \\
        --scored ${scored} \\
        --specificity ${spec_report} \\
        --clusters_species ${clusters_species} \\
        --manifest ${manifest} \\
        --all_cds ${all_cds} \\
        --rank species \\
        --cap ${params.marker_ani_cap} \\
        --max_copies ${params.marker_ani_max_copies} \\
        --kmer ${params.marker_ani_kmer} \\
        --ymax ${params.marker_ani_dist_max} \\
        --out marker_ani.json
    """

    stub:
    """
    printf '{"kmer":%s,"cap":%s,"max_copies":%s,"ymax":%s,"stages":["specific-200","post-crossmap","post-recovery"],"species":{}}' ${params.marker_ani_kmer} ${params.marker_ani_cap} ${params.marker_ani_max_copies} ${params.marker_ani_dist_max} > marker_ani.json
    """
}
