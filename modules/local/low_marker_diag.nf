// Diagnostics for species that ended below the marker threshold
// (params.low_marker_threshold). Two questions per flagged species:
//   1. can masking rescue the markers the cross-map guard dropped, or do they
//      cross-map whole-gene? (LOW_MARKER_MASKING -> assess via reps.ffn)
//   2. is the within-species ANI cleanly tighter than the between-species ANI in
//      the genus, or do they overlap (a re-merge signal)? (ANI_GAP -> skani)
// LOW_MARKER_REPORT folds both into one self-contained HTML page.

process LOW_MARKER_MASKING {
    tag "low_marker_masking"
    publishDir "${params.outdir}/low_marker", mode: 'copy'

    input:
    path manifest        // manifest.tsv
    path markers         // markers.specific.tsv (final, post-QC)
    path spec_report     // specificity_report.tsv
    path idmap           // query.map
    path hits            // crossmap.m8
    path marker_fasta    // markers.nuc.fasta (marker CDS / crossmap queries)
    path target          // reps.ffn (species-rep CDS; fast localisation target)

    output:
    path 'low_marker_species.tsv', emit: species
    path 'masking_summary.tsv',    emit: summary
    path 'masking_markers.tsv',    emit: markers
    path 'masking_coverage.txt',   emit: coverage

    script:
    """
    low_marker_masking.py \\
        --manifest ${manifest} \\
        --markers ${markers} \\
        --report ${spec_report} \\
        --idmap ${idmap} \\
        --crossmap ${hits} \\
        --marker_fasta ${marker_fasta} \\
        --target ${target} \\
        --threshold ${params.low_marker_threshold} \\
        --min_id ${params.specificity_min_id} \\
        --min_aln ${params.specificity_min_aln} \\
        --min_survivor ${params.masking_min_survivor} \\
        --outdir .
    """

    stub:
    """
    touch low_marker_species.tsv masking_summary.tsv masking_markers.tsv masking_coverage.txt
    """
}

process ANI_SKANI {
    tag "ani_skani"
    publishDir "${params.outdir}/low_marker", mode: 'copy'
    // Expensive all-vs-all triangle; hash inputs by path+size (not content+
    // timestamp) so -resume reuses it on shared/HPC filesystems.
    cache 'lenient'

    // skani runs in its own image (which has NO python) -- so this process does
    // ONLY the skani call; the python aggregation is a separate process below.
    conda "bioconda::skani=0.3.2"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/skani:0.3.2--h79ce301_0'
        : 'quay.io/biocontainers/skani:0.3.2--h79ce301_0'}"

    input:
    path genomes               // staged genome FASTAs for flagged genera

    output:
    path 'skani.sparse.tsv', emit: sparse

    script:
    // One all-vs-all triangle over every staged genome; the -s 80 screen drops
    // cross-genus pairs (<80% ANI) automatically, and the python step keeps only
    // same-genus pairs anyway. -E = sparse edge list (Ref Query ANI af_ref af_query).
    """
    find -L . -maxdepth 1 -type f \\
        \\( -name '*.fna.gz' -o -name '*.fa.gz' -o -name '*.fasta.gz' \\
           -o -name '*.fna' -o -name '*.fa' -o -name '*.fasta' \\) > genomes.txt
    skani triangle \\
        -l genomes.txt \\
        -o skani.sparse.tsv \\
        -E \\
        --min-af ${params.low_marker_min_af * 100} \\
        -s 80 \\
        -t ${task.cpus}
    """

    stub:
    """
    printf 'Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\n' > skani.sparse.tsv
    """
}

process ANI_GAP {
    tag "ani_gap"
    publishDir "${params.outdir}/low_marker", mode: 'copy'

    // Pure stdlib-python: parse the skani sparse table, split within/between,
    // compute the ANI gap. Runs in the default env (has python3), NOT the skani
    // container.

    input:
    path manifest              // manifest.tsv
    path low_marker_species    // low_marker_species.tsv (flagged column)
    path sparse                // skani.sparse.tsv from ANI_SKANI

    output:
    path 'ani_pairs.tsv',       emit: pairs
    path 'ani_gap_summary.tsv', emit: gap

    script:
    """
    ani_gap.py \\
        --manifest ${manifest} \\
        --low_marker_species ${low_marker_species} \\
        --sparse ${sparse} \\
        --ani ${params.low_marker_ani} \\
        --min_af ${params.low_marker_min_af} \\
        --outdir .
    """

    stub:
    """
    printf 'focal_species\\tother_species\\tkind\\tani\\talign_frac\\n' > ani_pairs.tsv
    printf 'species\\tgenus\\tn_within\\tmin_within\\tmedian_within\\tn_between\\tnearest_species\\tmax_between\\tgap\\toverlap\\tmerge_candidate\\n' > ani_gap_summary.tsv
    """
}

process MERGE_GAIN {
    tag "merge_gain"
    publishDir "${params.outdir}/low_marker", mode: 'copy'

    // Pure stdlib-python: for every flagged species, greedily merge its most
    // overlapping same-genus neighbours and recompute the marker set from
    // counts.tsv alone (merging is only relabelling -- no re-run), stopping once
    // the capped marker count reaches low_marker_threshold. Answers "would
    // merging lift this species to the requested marker count, and by how much?".

    input:
    path manifest              // manifest.tsv
    path counts                // counts.tsv
    path low_marker_species    // low_marker_species.tsv (flagged column)

    output:
    path 'merge_gain.tsv', emit: gain

    script:
    """
    merge_gain.py \\
        --manifest ${manifest} \\
        --counts ${counts} \\
        --low_marker_species ${low_marker_species} \\
        --threshold ${params.low_marker_threshold} \\
        --min_in ${params.min_in_prevalence} \\
        --max_out ${params.max_out_prevalence} \\
        --min_clade_size ${params.min_clade_size} \\
        --max_per_clade ${params.max_markers_per_clade} \\
        --outdir .
    """

    stub:
    """
    printf 'species\\tgenus\\tn_genomes\\tbaseline_markers\\tmerged_markers\\tdelta\\treached_threshold\\tn_species_added\\tmerged_clade\\ttrajectory\\n' > merge_gain.tsv
    """
}

process MASKING_REPORT {
    tag "masking_report"
    publishDir "${params.outdir}/low_marker", mode: 'copy'

    input:
    path masking_summary
    path masking_markers

    output:
    path 'masking_report.html', emit: report

    script:
    """
    masking_report.py \\
        --masking_summary ${masking_summary} \\
        --masking_markers ${masking_markers} \\
        --min_survivor ${params.masking_min_survivor} \\
        --out masking_report.html
    """

    stub:
    """
    touch masking_report.html
    """
}

process LOW_MARKER_REPORT {
    tag "low_marker_report"
    publishDir "${params.outdir}/low_marker", mode: 'copy'

    input:
    path low_marker_species
    path masking_summary
    path ani_pairs
    path ani_gap_summary
    path merge_gain

    output:
    path 'low_marker_report.html', emit: report

    script:
    """
    low_marker_report.py \\
        --low_marker_species ${low_marker_species} \\
        --masking_summary ${masking_summary} \\
        --ani_pairs ${ani_pairs} \\
        --ani_gap_summary ${ani_gap_summary} \\
        --merge_gain ${merge_gain} \\
        --threshold ${params.low_marker_threshold} \\
        --ani_threshold ${params.low_marker_ani} \\
        --out low_marker_report.html
    """

    stub:
    """
    touch low_marker_report.html
    """
}
