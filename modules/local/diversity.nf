process DIVERSITY {
    tag "diversity"
    publishDir "${params.outdir}/markers", mode: 'copy'

    input:
    path marker_fasta
    path markers        // scored marker table to augment

    output:
    path 'markers.diversity.tsv', emit: marker_table

    script:
    // Within-clade nucleotide divergence per marker via mash sketches of the
    // per-species-rep copies. Feeds #reps-per-marker and classifier threshold.
    """
    diversity.py \\
        --fasta ${marker_fasta} \\
        --markers ${markers} \\
        --kmer ${params.div_kmer} \\
        --threads ${task.cpus} \\
        --out markers.diversity.tsv
    """
}
