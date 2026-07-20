// Core-gene mask rescue for low-marker species. For every FLAGGED species, take
// ALL its core genes (in_prevalence >= min_in) with the uniqueness gate DROPPED,
// cross-map them against all-CDS, and keep the ones the cross-map guard can
// recover BY MASKING -- i.e. a long clean (non-cross-mapping) window survives even
// though the gene is shared. This targets markers that live in the locally
// divergent window of an otherwise-shared core gene, which the uniqueness filter
// discards (overall uniqueness and local maskability are orthogonal axes).
//
// Differs from the relax path in two ways: (1) it drops uniqueness instead of
// relaxing in-prevalence, and (2) CORE_GUARD runs specificity_guard.py WITH
// --target (mask recovery ON), whereas RELAX_GUARD is clean-only. Separate module
// so published intermediates land under markers/core/ and don't clash.

process SELECT_CORE {
    tag "select_core"
    publishDir "${params.outdir}/markers/core", mode: 'copy'

    input:
    path counts
    path clade_sizes
    path low_marker      // low_marker_species.tsv (flagged column)
    path emitted         // markers.emitted.tsv (exclude already-guarded)

    output:
    path 'core_candidates.tsv', emit: candidates

    script:
    """
    select_core.py \\
        --counts ${counts} \\
        --clade_sizes ${clade_sizes} \\
        --low_marker ${low_marker} \\
        --emitted ${emitted} \\
        --rank species \\
        --min_in ${params.min_in_prevalence} \\
        --min_clade_size ${params.min_clade_size} \\
        --score_out_exp ${params.score_out_exp} \\
        --max_per_clade ${params.core_rescue_max_per_clade} \\
        --out core_candidates.tsv
    """

    stub:
    "printf 'rank\\tclade\\tcluster\\tin_prevalence\\tout_prevalence\\tin_count\\tclade_size\\tscore\\n' > core_candidates.tsv"
}

process CORE_EMIT {
    tag "core_emit"
    publishDir "${params.outdir}/markers/core", mode: 'copy'

    input:
    path markers            // core_candidates.tsv
    path clusters_loose
    path clusters_species
    path reps_ffn
    path manifest

    output:
    path 'core.nuc.fasta',   emit: marker_fasta
    path 'core.emitted.tsv', emit: markers

    script:
    """
    emit_reps.py \\
        --markers ${markers} \\
        --clusters ${clusters_loose} \\
        --clusters_species ${clusters_species} \\
        --reps_ffn ${reps_ffn} \\
        --manifest ${manifest} \\
        --min_gene_len ${params.min_gene_len} \\
        --out_markers core.emitted.tsv \\
        --out core.nuc.fasta
    """

    stub:
    "touch core.nuc.fasta core.emitted.tsv"
}

process CORE_CROSSMAP {
    tag "core_crossmap"
    label 'process_high'
    publishDir "${params.outdir}/markers/core", mode: 'copy'
    cache 'lenient'

    conda "bioconda::mmseqs2=18.8cc5c"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"

    input:
    path marker_fasta   // core candidate CDS (queries)
    path target_cds     // all_cds.ffn (cross-map target universe)

    output:
    path 'core.m8',  emit: hits
    path 'core.map', emit: idmap

    script:
    def split_mem = task.memory ? "${(task.memory.toGiga() * 0.85) as long}G" : ''
    def split_arg = split_mem ? "--split-memory-limit ${split_mem}" : ''
    """
    awk 'BEGIN{n=0}
         /^>/ { n++; print "m"n"\\t"substr(\$0,2) > "core.map"; print ">m"n; next }
         { print }' ${marker_fasta} > core_query.fna

    # No flagged-species core candidates -> nothing to search; emit empty and stop.
    touch core.map
    if [ ! -s core_query.fna ]; then touch core.m8; exit 0; fi

    mmseqs easy-search \\
        core_query.fna \\
        ${target_cds} \\
        core.m8 \\
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
    "touch core.m8 core.map"
}

process CORE_GUARD {
    tag "core_guard"
    publishDir "${params.outdir}/markers/core", mode: 'copy'

    // Mask-recovery guard: specificity_guard.py WITH --target all_cds, so a shared
    // core gene is KEPT when a clean window >= recovery_min_clean survives with
    // <= recovery_max_masked_frac of the gene masked. This is the whole point of
    // the core-rescue path -- unlike RELAX_GUARD, cross-mapping is not fatal if a
    // long enough clean window remains.
    input:
    path hits          // core.m8
    path idmap         // core.map
    path markers       // core.emitted.tsv
    path marker_fasta  // core.nuc.fasta
    path manifest
    path target_cds    // all_cds.ffn (mask-recovery target universe)

    output:
    path 'core.recovered.tsv',   emit: markers
    path 'core.recovered.fasta', emit: marker_fasta
    path 'core.report.tsv',      emit: report

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
        --target ${target_cds} \\
        --recovery_min_clean ${params.mask_recovery_min_clean} \\
        --recovery_max_masked_frac ${params.mask_recovery_max_masked_frac} \\
        --out_markers core.recovered.tsv \\
        --out_fasta core.recovered.fasta \\
        --report core.report.tsv
    """

    stub:
    "touch core.recovered.tsv core.recovered.fasta core.report.tsv"
}

process SELECT_CORE_RESCUED {
    tag "select_core_rescued"
    publishDir "${params.outdir}/markers/core", mode: 'copy'

    input:
    path base_markers   // current final markers (base or post-relax)
    path extra          // core.recovered.tsv
    path extra_fasta    // core.recovered.fasta

    output:
    path 'core.rescued.markers.tsv', emit: markers
    path 'core.rescued.fasta',       emit: marker_fasta

    script:
    """
    select_relaxed.py \\
        --base_markers ${base_markers} \\
        --extra_clean ${extra} \\
        --extra_fasta ${extra_fasta} \\
        --rank species \\
        --threshold ${params.low_marker_threshold} \\
        --out_markers core.rescued.markers.tsv \\
        --out_fasta core.rescued.fasta
    """

    stub:
    "touch core.rescued.markers.tsv core.rescued.fasta"
}

process MERGE_CORE {
    tag "merge_core"
    publishDir "${params.outdir}/markers", mode: 'copy'

    input:
    path base_markers   // current final markers table
    path base_fasta     // current final marker FASTA
    path rescued_markers
    path rescued_fasta

    output:
    path 'markers.final.tsv',   emit: markers
    path 'markers.final.fasta', emit: marker_fasta

    script:
    // base already carries the header; append the rescued rows (skip their header).
    """
    cp ${base_markers} markers.final.tsv
    tail -n +2 ${rescued_markers} >> markers.final.tsv
    cat ${base_fasta} ${rescued_fasta} > markers.final.fasta
    """

    stub:
    "touch markers.final.tsv markers.final.fasta"
}
