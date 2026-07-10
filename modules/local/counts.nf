process COUNTS {
    tag "counts"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path clusters_loose    // rep <TAB> member  (loose protein clustering)
    path clusters_species  // rep <TAB> member  (species-tight protein clustering)
    path manifest          // idx / genome_id / is_rep / path / lineage...

    output:
    path 'counts.tsv',      emit: counts        // cluster rank clade in_count marker_total
    path 'clade_sizes.tsv', emit: clade_sizes   // rank clade size   (+ __TOTAL__ row)

    script:
    // External sort groups members by cluster so the aggregator holds only one
    // cluster's genome set in RAM at a time. This is the 800k-safe path.
    // Both clusterings are streamed: loose owns every rank above species, the
    // tight species clustering owns the species rank (see aggregate_counts.py).
    """
    LC_ALL=C sort -k1,1 --parallel=${task.cpus} -S 25% ${clusters_loose}   > clusters_loose.sorted.tsv
    LC_ALL=C sort -k1,1 --parallel=${task.cpus} -S 25% ${clusters_species} > clusters_species.sorted.tsv

    aggregate_counts.py \\
        --clusters clusters_loose.sorted.tsv \\
        --clusters_species clusters_species.sorted.tsv \\
        --manifest ${manifest} \\
        --counts counts.tsv \\
        --clade_sizes clade_sizes.tsv
    """
}
