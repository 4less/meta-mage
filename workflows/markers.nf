nextflow.enable.dsl = 2

// nf-core modules. linclust/createtsv run twice: a loose (family-and-above)
// clustering and a tight species clustering. See conf/modules.config for args.
include { PYRODIGAL         } from '../modules/nf-core/pyrodigal/main'
include { MMSEQS_CREATEDB   } from '../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_LINCLUST as MMSEQS_LINCLUST_LOOSE    } from '../modules/nf-core/mmseqs/linclust/main'
include { MMSEQS_LINCLUST as MMSEQS_LINCLUST_SPECIES  } from '../modules/nf-core/mmseqs/linclust/main'
include { MMSEQS_CREATETSV as MMSEQS_CREATETSV_LOOSE   } from '../modules/nf-core/mmseqs/createtsv/main'
include { MMSEQS_CREATETSV as MMSEQS_CREATETSV_SPECIES } from '../modules/nf-core/mmseqs/createtsv/main'

// local (bespoke) modules
include { BUILD_MANIFEST } from '../modules/local/build_manifest.nf'
include { REHEADER as REHEADER_FAA } from '../modules/local/reheader.nf'
include { REHEADER as REHEADER_FNA } from '../modules/local/reheader.nf'
include { COUNTS         } from '../modules/local/counts.nf'
include { SCORE; SCORE_RELAX } from '../modules/local/score.nf'
include { FILTER_RELAX; RELAX_EMIT; RELAX_CROSSMAP; RELAX_GUARD; SELECT_RELAXED; MERGE_RELAX } from '../modules/local/relax.nf'
include { EMIT_REPS      } from '../modules/local/emit_reps.nf'
include { CROSSMAP       } from '../modules/local/crossmap.nf'
include { SPECIFICITY    } from '../modules/local/specificity.nf'
include { NAKEDNESS_SEARCH; NAKEDNESS } from '../modules/local/nakedness.nf'
include { MARKER_ANI     } from '../modules/local/marker_ani.nf'
include { DIVERSITY      } from '../modules/local/diversity.nf'
include { EMIT_DB        } from '../modules/local/emit_db.nf'
include { REPORT; LEAKAGE_REPORT } from '../modules/local/report.nf'
include { LOW_MARKER_MASKING; ANI_SKANI; ANI_GAP; MERGE_GAIN; MASKING_REPORT; LOW_MARKER_REPORT } from '../modules/local/low_marker_diag.nf'

workflow MARKERS {

    // 1. Integer-indexed genome manifest with full lineage (single source of truth).
    //    Optional GTDB metadata attaches CheckM2 completeness for the weighted
    //    core-prevalence; NO_FILE placeholder keeps the input stable when unset.
    gtdb_meta = params.gtdb_metadata ? file(params.gtdb_metadata)
                                     : file("${projectDir}/assets/NO_FILE")
    BUILD_MANIFEST(file(params.genome_sheet), gtdb_meta)
    manifest = BUILD_MANIFEST.out.manifest

    // 2. Per-genome channel (the big SCATTER). meta carries the integer idx.
    genomes = manifest
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(
            [ id: row.genome_id, idx: row.idx, is_rep: row.is_rep == '1',
              genus: row.genus, species: row.species ],
            file(row.path)
        ) }

    // 3. Gene calling. nf-core pyrodigal emits BOTH proteins (.faa) and
    //    nucleotide CDS (.fna) in one pass.
    PYRODIGAL(genomes, 'gff')

    // 4. Make protein IDs globally unique (g<idx>_<orig>) so clustering can't
    //    confuse identically-named contigs across genomes.
    REHEADER_FAA(PYRODIGAL.out.faa)

    all_faa = REHEADER_FAA.out.fasta
        .map { meta, fa -> fa }
        .collectFile(name: 'all_proteins.faa', storeDir: "${params.outdir}/linclust")

    // 5. Two protein clusterings over the SAME protein DB -> rep<TAB>member TSVs
    //    (the vocabulary). Loose scaffolds family-and-above; tight (species)
    //    keeps sister-species / paralog genes apart for species specificity.
    MMSEQS_CREATEDB(all_faa.map { fa -> tuple([id: 'proteins'], fa) })

    MMSEQS_LINCLUST_LOOSE(MMSEQS_CREATEDB.out.db)
    MMSEQS_CREATETSV_LOOSE(
        MMSEQS_LINCLUST_LOOSE.out.db_cluster,
        MMSEQS_CREATEDB.out.db,
        MMSEQS_CREATEDB.out.db
    )
    clusters_loose = MMSEQS_CREATETSV_LOOSE.out.tsv.map { meta, tsv -> tsv }

    MMSEQS_LINCLUST_SPECIES(MMSEQS_CREATEDB.out.db)
    MMSEQS_CREATETSV_SPECIES(
        MMSEQS_LINCLUST_SPECIES.out.db_cluster,
        MMSEQS_CREATEDB.out.db,
        MMSEQS_CREATEDB.out.db
    )
    clusters_species = MMSEQS_CREATETSV_SPECIES.out.tsv.map { meta, tsv -> tsv }

    // 6. Presence/absence as per-clade COUNTS, streamed one cluster at a time.
    //    Species rank is counted from the tight clustering, higher ranks from
    //    the loose one; cluster ids are namespaced by source (L:/S:).
    COUNTS(clusters_loose, clusters_species, manifest)

    // 7. Specificity scoring across all ranks.
    SCORE(COUNTS.out.counts, COUNTS.out.clade_sizes)

    // 8. Nucleotide CDS, reheadered (g<idx>_<orig>) for EVERY genome. The marker
    //    payload is drawn from species reps only (reps_fna), but the specificity
    //    guard needs the full corpus (all_cds) so cross-mapping is tested against
    //    every genome, not just one representative per species.
    REHEADER_FNA(PYRODIGAL.out.fna)
    reheadered_fna = REHEADER_FNA.out.fasta

    reps_fna = reheadered_fna
        .filter { meta, fna -> meta.is_rep }
        .map { meta, fna -> fna }
        .collectFile(name: 'reps.ffn', storeDir: "${params.outdir}/emit")

    all_cds = reheadered_fna
        .map { meta, fna -> fna }
        .collectFile(name: 'all_cds.ffn', storeDir: "${params.outdir}/emit")

    // EMIT_REPS also drops markers whose every rep CDS is < min_gene_len and emits
    // the length-filtered marker table (markers.emitted.tsv); everything downstream
    // uses that so the QC set, the DB, and the report all agree.
    EMIT_REPS(SCORE.out.markers, clusters_loose, clusters_species, reps_fna, manifest)
    markers_emitted = EMIT_REPS.out.markers

    // 9. Nucleotide specificity guard (MetaPhlAn's final uniqueness check): align
    //    marker CDS against ALL genomes' CDS and drop markers that cross-map to
    //    an off-target clade at read-mapping identity.
    // Merge assessment for the final report; NO_FILE unless the low-marker stage runs.
    merge_gain_report = file("${projectDir}/assets/NO_FILE")
    // Nakedness of cross-map drops; NO_FILE unless the guard runs.
    nakedness_report = file("${projectDir}/assets/NO_FILE")
    if( params.specificity ) {
        CROSSMAP(EMIT_REPS.out.marker_fasta, all_cds)
        SPECIFICITY(CROSSMAP.out.hits, CROSSMAP.out.idmap, markers_emitted, EMIT_REPS.out.marker_fasta, manifest, all_cds)
        markers_final = SPECIFICITY.out.markers
        fasta_final   = SPECIFICITY.out.marker_fasta
        spec_report   = SPECIFICITY.out.report

        // Of the markers the guard dropped, how many are NAKED (no competing
        // off-target marker -> truly steal reads) vs CONTESTED (the off-target
        // clade markers the region too -> conservatively dropped). Marker-vs-
        // marker self search + classifier.
        NAKEDNESS_SEARCH(EMIT_REPS.out.marker_fasta)
        NAKEDNESS(NAKEDNESS_SEARCH.out.hits, NAKEDNESS_SEARCH.out.idmap,
                  SPECIFICITY.out.report)
        nakedness_report = NAKEDNESS.out.nakedness

        // Per-species cross-map leakage report (incoming/outgoing, dropdown focus).
        LEAKAGE_REPORT(SPECIFICITY.out.report)

        // 9b. Low-marker diagnostics: for species with < low_marker_threshold final
        //     markers, ask whether their dropped markers are rescuable by masking,
        //     and whether the species has a within/between ANI gap in its genus.
        if( params.low_marker_threshold > 0 ) {
            LOW_MARKER_MASKING(
                manifest,
                SPECIFICITY.out.markers,
                SPECIFICITY.out.report,
                CROSSMAP.out.idmap,
                CROSSMAP.out.hits,
                EMIT_REPS.out.marker_fasta,
                reps_fna
            )

            // Dedicated cross-map masking report (per-marker rescue verdicts).
            MASKING_REPORT(LOW_MARKER_MASKING.out.summary,
                           LOW_MARKER_MASKING.out.markers)

            // Genera that contain at least one flagged species.
            flagged_genera = LOW_MARKER_MASKING.out.species
                .splitCsv(header: true, sep: '\t')
                .filter { it.flagged == 'yes' }
                .map { it.genus }
                .unique()
                .collect()

            // Stage only the genomes of those genera for skani. toSortedList (not
            // collect) so the staged file order is deterministic across runs --
            // otherwise the task hash changes every time and -resume re-runs skani.
            ani_genomes = genomes
                .combine(flagged_genera)
                .filter { meta, fa, genera -> genera.contains(meta.genus) }
                .map { meta, fa, genera -> fa }
                .toSortedList()

            ANI_SKANI(ani_genomes)
            ANI_GAP(manifest, LOW_MARKER_MASKING.out.species, ANI_SKANI.out.sparse)

            // Does merging a flagged species with its most overlapping same-genus
            // neighbours lift it to low_marker_threshold markers? Recomputed from
            // counts.tsv alone (no re-run) -- the metric for whether a merge pays off.
            MERGE_GAIN(manifest, COUNTS.out.counts, LOW_MARKER_MASKING.out.species,
                       SPECIFICITY.out.report, SPECIFICITY.out.mask_intervals)
            merge_gain_report = MERGE_GAIN.out.gain

            LOW_MARKER_REPORT(
                LOW_MARKER_MASKING.out.species,
                LOW_MARKER_MASKING.out.summary,
                ANI_GAP.out.pairs,
                ANI_GAP.out.gap,
                MERGE_GAIN.out.gain
            )

            // 9c. Adaptive in-prevalence relaxation. For species still below the
            //     threshold after mask recovery AND without a threshold-reaching
            //     merge (FILTER_RELAX picks them), lower the in-prevalence floor,
            //     guard the extra candidates against all-CDS, keep ONLY cross-map-
            //     free markers, and add higher-prevalence tiers first until the
            //     threshold is met. Rescued markers fold into the final DB. If
            //     nothing is flagged the branch runs empty and the output equals
            //     the base set.
            if( params.relax_in_prevalence ) {
                SCORE_RELAX(COUNTS.out.counts, COUNTS.out.clade_sizes)
                FILTER_RELAX(SCORE_RELAX.out.markers,
                             LOW_MARKER_MASKING.out.species, MERGE_GAIN.out.gain)
                RELAX_EMIT(FILTER_RELAX.out.candidates,
                           clusters_loose, clusters_species, reps_fna, manifest)
                RELAX_CROSSMAP(RELAX_EMIT.out.marker_fasta, all_cds)
                RELAX_GUARD(RELAX_CROSSMAP.out.hits, RELAX_CROSSMAP.out.idmap,
                            RELAX_EMIT.out.markers, RELAX_EMIT.out.marker_fasta,
                            manifest)
                SELECT_RELAXED(SPECIFICITY.out.markers,
                               RELAX_GUARD.out.markers, RELAX_GUARD.out.marker_fasta)
                MERGE_RELAX(SPECIFICITY.out.markers, SPECIFICITY.out.marker_fasta,
                            SELECT_RELAXED.out.markers, SELECT_RELAXED.out.marker_fasta)
                markers_final = MERGE_RELAX.out.markers
                fasta_final   = MERGE_RELAX.out.marker_fasta
            }
        }
    } else {
        markers_final = markers_emitted
        fasta_final   = EMIT_REPS.out.marker_fasta
        // No cross-map QC ran: stage a placeholder so REPORT's input is stable.
        spec_report   = file("${projectDir}/assets/NO_FILE")
    }

    // 10. Within-clade nucleotide divergence per marker (augments the marker table).
    DIVERSITY(fasta_final, markers_final)

    // 11. Final classifier DB: FASTA + annotated marker table.
    EMIT_DB(DIVERSITY.out.marker_table, fasta_final)

    // 11b. Per-species pairwise marker ANI at each filtering stage (report
    //      boxplots): for each top-scoring marker, pairwise distance among its
    //      copies across the species' genomes. Uses the pre-QC emitted marker set
    //      so the specific-N stage is complete; the guard verdicts split it into
    //      post-crossmap / -recovery. Copies come from the all-CDS set via the
    //      species clusters.
    MARKER_ANI(markers_emitted, spec_report, clusters_species, manifest, all_cds)

    // 12. Human-readable HTML report: the pangenome->core->marker funnel per
    //     species, why genes were removed at each QC stage, and which species
    //     were dropped (per-genus). Reads the pipeline's own tables so the
    //     numbers match the emitted DB.
    REPORT(
        manifest,
        BUILD_MANIFEST.out.dropped,
        COUNTS.out.counts,
        COUNTS.out.clade_sizes,
        markers_emitted,
        spec_report,
        merge_gain_report,
        nakedness_report,
        MARKER_ANI.out.ani
    )
}
