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
    path cross_hits     // crossmap.m8 (marker-vs-all_cds hits)
    path cross_map      // query.map (m<N> -> rank|clade|cluster|genome)
    path manifest       // idx -> lineage (off-target genome clade lookup)
    path spec_report    // specificity_report.tsv (pass column marks guard-dropped)

    output:
    path 'nakedness.tsv', emit: nakedness

    script:
    // Of the markers the guard dropped, which are NAKED (off-target region maps
    // BETTER to this marker than to the off-target clade's own competing marker ->
    // truly steals reads) vs CONTESTED (a competitor matches >= as well -> reads
    // stay home, drop was conservative)? Compares crossmap id_off vs self-search
    // id_comp per off-target clade.
    """
    nakedness.py \\
        --self_hits ${self_hits} \\
        --selfmap ${self_map} \\
        --cross_hits ${cross_hits} \\
        --crossmap ${cross_map} \\
        --manifest ${manifest} \\
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


// Competitive rescue: re-admit the guard-dropped markers NAKEDNESS labelled
// CONTESTED (every off-target clade has an >= competing marker). Their rows and
// sequences are pulled from the PRE-guard emitted set (the guard removed them from
// markers.specific.*), then MERGE_NAKEDNESS unions them onto the guard survivors.
process NAKEDNESS_SELECT {
    tag "nakedness_select"
    publishDir "${params.outdir}/markers", mode: 'copy'

    input:
    path nakedness      // nakedness.tsv (verdict column)
    path emitted        // markers.emitted.tsv (pre-guard marker table)
    path emitted_fasta  // markers.nuc.fasta (pre-guard sequences)

    output:
    path 'contested.markers.tsv', emit: markers
    path 'contested.fasta',       emit: marker_fasta

    script:
    def rank_arg = params.competitive_rescue_rank ? "--rank ${params.competitive_rescue_rank}" : ''
    """
    select_contested.py \\
        --nakedness ${nakedness} \\
        --markers ${emitted} \\
        --fasta ${emitted_fasta} \\
        ${rank_arg} \\
        --out_markers contested.markers.tsv \\
        --out_fasta contested.fasta
    """

    stub:
    "touch contested.markers.tsv contested.fasta"
}


process MERGE_NAKEDNESS {
    tag "merge_nakedness"
    publishDir "${params.outdir}/markers", mode: 'copy'

    input:
    path base_markers    // markers.specific.tsv (guard survivors)
    path base_fasta      // markers.specific.fasta
    path rescued_markers // contested.markers.tsv
    path rescued_fasta   // contested.fasta

    output:
    path 'markers.compet.tsv',   emit: markers
    path 'markers.compet.fasta', emit: marker_fasta

    script:
    // base already has the header; append the rescued rows (skip their header).
    """
    cp ${base_markers} markers.compet.tsv
    tail -n +2 ${rescued_markers} >> markers.compet.tsv
    cat ${base_fasta} ${rescued_fasta} > markers.compet.fasta
    """

    stub:
    "touch markers.compet.tsv markers.compet.fasta"
}
