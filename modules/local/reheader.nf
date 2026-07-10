process REHEADER {
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("g${meta.idx}.reheader.fasta"), emit: fasta

    script:
    // Prefix every header with the integer genome id: g<idx>_<original-token>.
    // Keeps ids globally unique AND compact; faa/fna stay in correspondence
    // because we prefix (not renumber) the original pyrodigal id.
    """
    zcat -f ${fasta} \\
        | awk -v idx=${meta.idx} '/^>/{ sub(/^>/,""); print ">g"idx"_"\$1; next } { print }' \\
        > g${meta.idx}.reheader.fasta
    """
}
