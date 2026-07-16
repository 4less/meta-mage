process PREVALENCE_CURVES {
    tag "prevalence_curves"
    publishDir "${params.outdir}/report", mode: 'copy'

    // stdlib-python only; base env like REPORT. Standalone HTML with a species
    // dropdown: per-species gene-prevalence rank curve over the pangenome.
    input:
    path counts        // counts.tsv
    path clade_sizes   // clade_sizes.tsv

    output:
    path 'prevalence_curves.html', emit: report

    script:
    """
    prevalence_curves.py \\
        --counts ${counts} \\
        --clade_sizes ${clade_sizes} \\
        --rank species \\
        --core_prev ${params.core_prevalence} \\
        --relax_floor ${params.relax_in_floor} \\
        --out prevalence_curves.html
    """

    stub:
    "touch prevalence_curves.html"
}
