process SPECIFICITY {
    tag "specificity"
    publishDir "${params.outdir}/markers", mode: 'copy'

    input:
    path hits          // crossmap.m8: query target pident alnlen
    path idmap         // query.map: m<N> <TAB> original header
    path markers       // scored marker table
    path marker_fasta  // emitted marker CDS
    path manifest
    path target_cds    // all_cds.ffn: off-target CDS universe for mask recovery

    output:
    path 'markers.specific.tsv',   emit: markers
    path 'markers.specific.fasta', emit: marker_fasta
    path 'specificity_report.tsv', emit: report
    path 'mask_intervals.tsv',     emit: mask_intervals

    script:
    // A cross-mapping marker is RECOVERED (kept whole) if masking the shared
    // stretches still leaves a clean window >= mask_recovery_min_clean bp with
    // <= mask_recovery_max_masked_frac of the gene masked. --target off = legacy
    // binary guard. mask_intervals.tsv records per-off-target-clade masked spans
    // so the merge probe can re-mask after a hypothetical merge.
    def recover = params.mask_recovery ? "--target ${target_cds} --recovery_min_clean ${params.mask_recovery_min_clean} --recovery_max_masked_frac ${params.mask_recovery_max_masked_frac}" : ''
    """
    specificity_guard.py \\
        --hits ${hits} \\
        --idmap ${idmap} \\
        --markers ${markers} \\
        --fasta ${marker_fasta} \\
        --manifest ${manifest} \\
        --min_id ${params.specificity_min_id} \\
        --min_aln ${params.specificity_min_aln} \\
        ${recover} \\
        --out_markers markers.specific.tsv \\
        --out_fasta markers.specific.fasta \\
        --report specificity_report.tsv \\
        --out_mask_intervals mask_intervals.tsv
    """
}
