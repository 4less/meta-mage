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
        --score_out_exp ${params.score_out_exp} \\
        --out markers.tsv
    """
}

// Same scoring, but at the relaxed in-prevalence floor and a larger per-clade cap,
// so the lower-prevalence tiers survive for the adaptive rescue (FILTER_RELAX then
// picks the flagged species' sub-0.8 candidates).
process SCORE_RELAX {
    tag "score_relax"

    input:
    path counts
    path clade_sizes

    output:
    path 'markers_relaxed.tsv', emit: markers

    script:
    """
    score_markers.py \\
        --counts ${counts} \\
        --clade_sizes ${clade_sizes} \\
        --min_in ${params.relax_in_floor} \\
        --max_out ${params.max_out_prevalence} \\
        --min_clade_size ${params.min_clade_size} \\
        --max_per_clade ${params.relax_max_per_clade} \\
        --score_out_exp ${params.score_out_exp} \\
        --out markers_relaxed.tsv
    """
}
