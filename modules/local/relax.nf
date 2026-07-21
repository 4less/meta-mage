// Adaptive in-prevalence relaxation: rescue markers for species still below the
// threshold after mask recovery and merge. Lower the in-prevalence floor, guard
// the extra candidates against all-CDS, keep ONLY cross-map-free markers, and add
// higher-prevalence tiers first until the threshold is met. Separate modules (not
// reused CROSSMAP/EMIT_REPS) so their published outputs don't clash with the main
// run's; intermediates go under markers/relax/.

process FILTER_RELAX {
    tag "filter_relax"
    publishDir "${params.outdir}/markers/relax", mode: 'copy'

    input:
    path relaxed        // markers_relaxed.tsv (scored at the floor)
    path low_marker     // low_marker_species.tsv
    path merge_gain     // merge_gain.tsv or NO_FILE
    path emitted        // markers.emitted.tsv (what the main run actually guarded)

    output:
    path 'relax_candidates.tsv', emit: candidates

    script:
    """
    filter_relax.py \\
        --relaxed ${relaxed} \\
        --low_marker ${low_marker} \\
        --merge_gain ${merge_gain} \\
        --emitted ${emitted} \\
        --rank species \\
        --out relax_candidates.tsv
    """

    stub:
    "printf 'rank\\tclade\\tcluster\\tin_prevalence\\tout_prevalence\\tin_count\\tclade_size\\tscore\\n' > relax_candidates.tsv"
}

process RELAX_EMIT {
    tag "relax_emit"
    publishDir "${params.outdir}/markers/relax", mode: 'copy'

    input:
    path markers           // relax_candidates.tsv
    path clusters_loose
    path clusters_species
    path all_cds
    path manifest

    output:
    path 'relax.nuc.fasta',   emit: marker_fasta
    path 'relax.emitted.tsv', emit: markers

    script:
    """
    emit_reps.py \\
        --markers ${markers} \\
        --clusters ${clusters_loose} \\
        --clusters_species ${clusters_species} \\
        --all_cds ${all_cds} \\
        --manifest ${manifest} \\
        --min_gene_len ${params.min_gene_len} \\
        --out_markers relax.emitted.tsv \\
        --out relax.nuc.fasta
    """

    stub:
    "touch relax.nuc.fasta relax.emitted.tsv"
}

process RELAX_CROSSMAP {
    tag "relax_crossmap"
    label 'process_high'
    publishDir "${params.outdir}/markers/relax", mode: 'copy'
    cache 'lenient'

    conda "bioconda::mmseqs2=18.8cc5c"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    path marker_fasta   // relaxed candidate CDS (queries)
    path target_cds     // all_cds.ffn (cross-map target universe)

    output:
    path 'relax.m8',  emit: hits
    path 'relax.map', emit: idmap

    script:
    def split_mem = task.memory ? "${(task.memory.toGiga() * 0.85) as long}G" : ''
    def split_arg = split_mem ? "--split-memory-limit ${split_mem}" : ''
    """
    awk 'BEGIN{n=0}
         /^>/ { n++; print "m"n"\\t"substr(\$0,2) > "relax.map"; print ">m"n; next }
         { print }' ${marker_fasta} > relax_query.fna

    # No flagged-species candidates -> nothing to search; emit empty and stop.
    touch relax.map
    if [ ! -s relax_query.fna ]; then touch relax.m8; exit 0; fi

    mmseqs easy-search \\
        relax_query.fna \\
        ${target_cds} \\
        relax.m8 \\
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
    "touch relax.m8 relax.map"
}

process RELAX_GUARD {
    tag "relax_guard"
    publishDir "${params.outdir}/markers/relax", mode: 'copy'

    // Clean-only guard: specificity_guard.py WITHOUT --target, so mask recovery is
    // off and any marker with an off-target hit is dropped. Cross-mapping is fully
    // eliminated from the rescued set.
    input:
    path hits          // relax.m8
    path idmap         // relax.map
    path markers       // relax.emitted.tsv
    path marker_fasta  // relax.nuc.fasta
    path manifest

    output:
    path 'relax.clean.tsv',   emit: markers
    path 'relax.clean.fasta', emit: marker_fasta

    script:
    """
    specificity_guard.py \\
        --hits ${hits} \\
        --idmap ${idmap} \\
        --markers ${markers} \\
        --fasta ${marker_fasta} \\
        --manifest ${manifest} \\
        --min_id ${params.specificity_min_id} \\
        --min_aln ${params.specificity_min_aln} \\
        --out_markers relax.clean.tsv \\
        --out_fasta relax.clean.fasta \\
        --report relax.clean.report.tsv
    """

    stub:
    "touch relax.clean.tsv relax.clean.fasta"
}

process SELECT_RELAXED {
    tag "select_relaxed"
    publishDir "${params.outdir}/markers/relax", mode: 'copy'

    input:
    path base_markers   // markers.specific.tsv (0.8 final)
    path extra_clean    // relax.clean.tsv
    path extra_fasta    // relax.clean.fasta

    output:
    path 'rescued.markers.tsv', emit: markers
    path 'rescued.fasta',       emit: marker_fasta

    script:
    """
    select_relaxed.py \\
        --base_markers ${base_markers} \\
        --extra_clean ${extra_clean} \\
        --extra_fasta ${extra_fasta} \\
        --rank species \\
        --threshold ${params.low_marker_threshold} \\
        --out_markers rescued.markers.tsv \\
        --out_fasta rescued.fasta
    """

    stub:
    "touch rescued.markers.tsv rescued.fasta"
}

process MERGE_RELAX {
    tag "merge_relax"
    publishDir "${params.outdir}/markers", mode: 'copy'

    input:
    path base_markers   // markers.specific.tsv
    path base_fasta     // markers.specific.fasta
    path rescued_markers
    path rescued_fasta

    output:
    path 'markers.final.tsv',   emit: markers
    path 'markers.final.fasta', emit: marker_fasta

    script:
    // base already has the header; append the rescued rows (skip their header).
    """
    cp ${base_markers} markers.final.tsv
    tail -n +2 ${rescued_markers} >> markers.final.tsv
    cat ${base_fasta} ${rescued_fasta} > markers.final.fasta
    """

    stub:
    "touch markers.final.tsv markers.final.fasta"
}
