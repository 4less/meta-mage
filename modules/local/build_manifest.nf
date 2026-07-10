process BUILD_MANIFEST {
    tag "manifest"
    publishDir "${params.outdir}/manifest", mode: 'copy'

    input:
    path genome_sheet   // accession, gtdb_lineage, path

    output:
    path 'manifest.tsv', emit: manifest

    script:
    """
    build_manifest.py \\
        --genome_sheet ${genome_sheet} \\
        --min_genomes_per_species ${params.min_genomes_per_species} \\
        --out manifest.tsv
    """
}
