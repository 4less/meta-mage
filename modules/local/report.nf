process REPORT {
    tag "report"
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path manifest          // manifest.tsv (kept genomes + lineage)
    path dropped           // dropped_species.tsv (prefilter drops)
    path counts            // counts.tsv (per-clade presence)
    path clade_sizes       // clade_sizes.tsv
    path scored            // markers.tsv (SCORE selected, post-cap)
    path spec_report       // specificity_report.tsv, or the NO_FILE placeholder

    output:
    path 'report.html', emit: report

    script:
    // stdlib-python only. The specificity report is optional: when the guard is
    // off, a NO_FILE placeholder is staged and the report omits the QC stage.
    """
    build_report.py \\
        --manifest ${manifest} \\
        --dropped ${dropped} \\
        --counts ${counts} \\
        --clade_sizes ${clade_sizes} \\
        --markers ${scored} \\
        --specificity ${spec_report} \\
        --min_in ${params.min_in_prevalence} \\
        --max_out ${params.max_out_prevalence} \\
        --min_clade_size ${params.min_clade_size} \\
        --max_per_clade ${params.max_markers_per_clade} \\
        --core_prevalence ${params.core_prevalence} \\
        --min_genomes_per_species ${params.min_genomes_per_species} \\
        --out report.html
    """

    stub:
    """
    touch report.html
    """
}
