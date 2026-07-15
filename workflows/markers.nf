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
include { SCORE          } from '../modules/local/score.nf'
include { EMIT_REPS      } from '../modules/local/emit_reps.nf'
include { CROSSMAP       } from '../modules/local/crossmap.nf'
include { SPECIFICITY    } from '../modules/local/specificity.nf'
include { DIVERSITY      } from '../modules/local/diversity.nf'
include { EMIT_DB        } from '../modules/local/emit_db.nf'
include { REPORT         } from '../modules/local/report.nf'
include { LOW_MARKER_MASKING; ANI_SKANI; ANI_GAP; MASKING_REPORT; LOW_MARKER_REPORT } from '../modules/local/low_marker_diag.nf'

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
    if( params.specificity ) {
        CROSSMAP(EMIT_REPS.out.marker_fasta, all_cds)
        SPECIFICITY(CROSSMAP.out.hits, CROSSMAP.out.idmap, markers_emitted, EMIT_REPS.out.marker_fasta, manifest)
        markers_final = SPECIFICITY.out.markers
        fasta_final   = SPECIFICITY.out.marker_fasta
        spec_report   = SPECIFICITY.out.report

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

            // Stage only the genomes of those genera for skani.
            ani_genomes = genomes
                .combine(flagged_genera)
                .filter { meta, fa, genera -> genera.contains(meta.genus) }
                .map { meta, fa, genera -> fa }
                .collect()

            ANI_SKANI(ani_genomes)
            ANI_GAP(manifest, LOW_MARKER_MASKING.out.species, ANI_SKANI.out.sparse)
            LOW_MARKER_REPORT(
                LOW_MARKER_MASKING.out.species,
                LOW_MARKER_MASKING.out.summary,
                ANI_GAP.out.pairs,
                ANI_GAP.out.gap
            )
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
        spec_report
    )
}
