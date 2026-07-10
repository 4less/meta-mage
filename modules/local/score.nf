process SCORE {
    tag "score"

    input:
    path counts        // cluster rank clade in_count marker_total
    path clade_sizes   // rank clade size (+ __TOTAL__)

    output:
    path 'markers.tsv', emit: markers

    script:
    """
    score_markers.py \\
        --counts ${counts} \\
        --clade_sizes ${clade_sizes} \\
        --min_in ${params.min_in_prevalence} \\
        --max_out ${params.max_out_prevalence} \\
        --min_clade_size ${params.min_clade_size} \\
        --max_per_clade ${params.max_markers_per_clade} \\
        --out markers.tsv
    """
}
