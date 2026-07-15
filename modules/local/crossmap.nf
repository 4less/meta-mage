process CROSSMAP {
    tag "crossmap"
    label 'process_high'
    publishDir "${params.outdir}/markers", mode: 'copy'
    // Expensive; on shared/HPC filesystems hash inputs by path+size, not content+
    // timestamp, so -resume doesn't spuriously re-run it after unrelated changes.
    cache 'lenient'

    // Same mmseqs2 image the nf-core modules use, so no extra provisioning.
    conda "bioconda::mmseqs2=18.8cc5c"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    path marker_fasta   // emitted marker CDS (queries)
    path target_cds     // CDS of ALL genomes (cross-map target universe)

    output:
    path 'crossmap.m8', emit: hits
    path 'query.map',   emit: idmap

    script:
    // GTDB clade names contain spaces (e.g. "s__Bacteroides ovatus"), and mmseqs
    // truncates sequence ids at whitespace -- which would corrupt the marker key
    // encoded in the header. Re-key each query to a compact whitespace-free id
    // (m<N>) with a side map back to the original rank|clade|cluster|genome_id.
    //
    // Nucleotide-vs-nucleotide (--search-type 3). -c 0 keeps local partial
    // alignments so a read-length window of cross-identity is still reported;
    // the alignment-length cut is applied downstream, not here.
    //
    // NOTE (scale): the target is now every genome, so in-clade hits (which we
    // discard) can crowd the per-query result set. --max-seqs must stay >= the
    // largest in-clade genome count for the guarded ranks, or genuine off-target
    // hits get pushed out and non-specific markers slip through. At GTDB scale
    // prebuild the target DB once (mmseqs createdb) instead of per-run easy-search.
    //
    // Bound mmseqs RAM to the memory SLURM actually granted this task: it splits
    // the target DB into chunks that fit the budget instead of loading it whole,
    // so the step can't OOM-kill the node. We hand it ~85% of task.memory and
    // keep the rest as headroom for the prefilter result set + process overhead.
    // Scales automatically with retries (task.memory grows per attempt).
    def split_mem = task.memory ? "${(task.memory.toGiga() * 0.85) as long}G" : ''
    def split_arg = split_mem ? "--split-memory-limit ${split_mem}" : ''
    """
    awk 'BEGIN{n=0}
         /^>/ { n++; print "m"n"\\t"substr(\$0,2) > "query.map"; print ">m"n; next }
         { print }' ${marker_fasta} > query.fna

    mmseqs easy-search \\
        query.fna \\
        ${target_cds} \\
        crossmap.m8 \\
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
    touch crossmap.m8 query.map
    """
}
