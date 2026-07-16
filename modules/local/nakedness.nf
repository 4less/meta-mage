// Nakedness of cross-map drops, in two steps mirroring CROSSMAP -> SPECIFICITY:
// the marker-vs-marker mmseqs self search runs in the mmseqs container, then the
// (stdlib-python) classifier runs in the base env -- the mmseqs image has no
// python, so the two must be separate processes.

process NAKEDNESS_SEARCH {
    tag "nakedness_search"
    publishDir "${params.outdir}/markers", mode: 'copy'
    cache 'lenient'

    // Same mmseqs2 image as CROSSMAP.
    conda "bioconda::mmseqs2=18.8cc5c"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    path marker_fasta   // emitted marker CDS (rank|clade|cluster|genome headers)

    output:
    path 'self.m8',  emit: hits
    path 'self.map', emit: idmap

    script:
    // Rekey the marker fasta to whitespace-free m<N> ids (clades contain spaces;
    // mmseqs truncates at whitespace) and search it against itself. Same identity
    // /window as the guard so the two agree.
    def split_mem = task.memory ? "${(task.memory.toGiga() * 0.85) as long}G" : ''
    def split_arg = split_mem ? "--split-memory-limit ${split_mem}" : ''
    """
    awk 'BEGIN{n=0}
         /^>/ { n++; print "m"n"\\t"substr(\$0,2) > "self.map"; print ">m"n; next }
         { print }' ${marker_fasta} > markers.fna

    mmseqs easy-search \\
        markers.fna \\
        markers.fna \\
        self.m8 \\
        tmp \\
        --search-type 3 \\
        --min-seq-id ${params.specificity_min_id} \\
        -c 0 \\
        --max-seqs ${params.specificity_max_seqs} \\
        ${split_arg} \\
        --format-output query,target,pident,alnlen \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch self.m8 self.map
    """
}

process NAKEDNESS {
    tag "nakedness"
    publishDir "${params.outdir}/markers", mode: 'copy'

    // stdlib-python only; runs in the base env like SPECIFICITY (the mmseqs image
    // has no python), so no container is declared.

    input:
    path self_hits      // self.m8 (marker-vs-marker mmseqs hits)
    path self_map       // self.map (m<N> -> rank|clade|cluster|genome)
    path spec_report    // specificity_report.tsv (needs offtarget_clades column)

    output:
    path 'nakedness.tsv', emit: nakedness

    script:
    // Of the markers the guard dropped, how many are NAKED (no competing off-
    // target marker -> truly steal reads) vs CONTESTED (the off-target clade
    // markers the region too -> conservatively dropped)?
    """
    nakedness.py \\
        --self_hits ${self_hits} \\
        --selfmap ${self_map} \\
        --specificity ${spec_report} \\
        --min_id ${params.specificity_min_id} \\
        --min_aln ${params.specificity_min_aln} \\
        --out nakedness.tsv
    """

    stub:
    """
    printf 'rank\\tclade\\tcluster\\tn_offtarget_clades\\tn_contested\\tn_naked\\tverdict\\tnaked_clades\\n' > nakedness.tsv
    """
}
