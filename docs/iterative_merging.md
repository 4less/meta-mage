# Iterative phylogeny-guided merging — roadmap

Recovering usable markers for species that share (nearly) all their gene content
with a tight relative complex (the motivating case: `s__Bacteroides ovatus`), by
merging clades **along the phylogeny** and **re-evaluating cross-mapping from
scratch** at each step — not by relabelling recorded results.

Status legend: ✅ done · 🚧 in progress · ⬜ not started.

---

## 1. Why

A species in a dense complex can end with **zero** specific markers even though it
is a valid GTDB species. For `B. ovatus` (v2 run):

- **82 / 82** candidate markers cross-map off-target (0 clean).
- **49** leak at **100 % nucleotide identity** (identical copies — no signal),
  31 at 99–<100 %, only **2** at 95–<99 % (a maskable window exists).
- Leakage is **poly­phyletic**: ovatus leaks to ~23 species / from ~20, spread
  across the genus (all within `g__Bacteroides`; 0 cross-genus). It does **not**
  follow the bac120 phylogeny.

The leakage signal and the phylogeny disagree, so a single clade-merge cannot
match the leakage graph, and the greedy core-overlap merge in
[`bin/merge_gain.py`](../bin/merge_gain.py) produces a taxonomy-breaking,
**polyphyletic** "+3 species" result (see the low-marker report).

### The trap we measured

[`bin/phylo_merge.py`](../bin/phylo_merge.py) walks the bac120 tree upward from a
focus species and, at each ancestral node, treats the whole monophyletic clade as
one taxon. For ovatus, leakage internalisation climbs cleanly with clade size:

| clade size | out-leak removed | in-leak removed | markers |
|---:|---:|---:|---:|
| 4  | 33 % | 42 % | 0 |
| 5  | 82 % | 76 % | 1 |
| 6  | 96 % | 82 % | 0 |
| 13 | 99 % | 100 % | 0 |
| 43 (whole genus) | 100 % | 100 % | 0 |

**Leakage disappears but markers never return.** This is only *partly* real. The
probe recomputes prevalence and re-masks residual off-targets, but it is a
**conservative lower bound**: it can only judge genes that were *emitted and
cross-mapped in the original run*. Genes that become clade-specific **only after**
the merge were excluded up front by `out_prev ≤ 0.02`, so they were never emitted,
never cross-mapped, and have no mask intervals — the probe drops all of them,
and they are the **biggest potential recovery source**.

**Conclusion:** to know whether merging recovers markers, we must re-run
emit → cross-map → guard/mask for the merged clade against the *out-of-clade*
universe. That is what this roadmap builds.

---

## 2. Method

For a candidate merged clade `C` (a monophyletic set of run species):

1. **Pick a representative genome** — the **medoid**: the genome with maximum mean
   ANI to all other genomes in `C` (minimises distance to the cluster). Source:
   `low_marker/skani.sparse.tsv` (already emitted; genome-level ANI incl. the
   ovatus complex). Fallback if a candidate cluster is absent from the medoid: use
   any clade member's CDS so the candidate set stays complete.
2. **Recompute the prevalence candidate pool** for `C` from `counts.tsv`
   (`in_prev ≥ min_in`, `out_prev ≤ max_out` over the merged clade). This is where
   previously-too-widespread shared-core genes become eligible.
3. **Re-emit** those candidates' CDS from the medoid; relabel all members of `C`
   as one clade.
4. **Fresh cross-map** the re-emitted CDS against `all_cds.ffn` **minus the
   in-clade genomes** — real alignments, not recycled `crossmap.m8` rows.
5. **Guard + re-mask** on the fresh cross-map: drop whole-gene off-target hits,
   mask partial ones, keep a marker iff a clean window `≥ recovery_min_clean`
   survives with `≤ recovery_max_masked_frac` masked. Because most former leakage
   sources are now in-clade, masking is recomputed against far fewer off-targets
   and can recover genes that previously masked to death.

Steps 1–2 are local and cheap. Steps 3–5 need `mmseqs` + genome CDS and belong on
the cluster as a Nextflow process — a near-clone of the relax subworkflow
([`modules/local/relax.nf`](../modules/local/relax.nf):
`RELAX_EMIT → RELAX_CROSSMAP → RELAX_GUARD`) with a medoid rep and clade
relabelling.

### Why this can genuinely recover ovatus (not just be "more correct")

The species ovatus is 100 %-identical to — `sp900755095`, `xylanisolvens` — are
exactly the ones a modest complex-merge pulls **in-clade**. Once they stop
counting as off-target, the residual out-of-clade hits may sit at 96–99 % with
maskable windows. The floor of 0 was partly an artifact of never testing the
shared-core genes.

---

## 3. What already exists ✅

| Piece | File | Role |
|---|---|---|
| Bacteroides subclade | `data/phylogeny/bacteroides_subclade.tree` | 250-tip clade from bac120_r232 |
| accession→species | `data/phylogeny/bacteroides_species.tsv` | tip labels |
| Leakage extractor | [`bin/leakage_edges.py`](../bin/leakage_edges.py) | crossmap.m8 → directed species edges |
| Leakage edges (v2) | `data/phylogeny/bacteroides_leakage_edges.tsv` | 432 directed edges |
| Leakage viz | [`bin/leakage_tree.py`](../bin/leakage_tree.py) | radial tree + in/out arcs, per-focus |
| Merge probe (lower bound) | [`bin/phylo_merge.py`](../bin/phylo_merge.py) | phylogeny walk, leakage% + markers per node |
| Greedy merge (legacy) | [`bin/merge_gain.py`](../bin/merge_gain.py) | core-overlap greedy; polyphyletic |
| Re-emit template | [`modules/local/relax.nf`](../modules/local/relax.nf) | EMIT→CROSSMAP→GUARD to clone |

Inputs on M3 (`output/v2/`): `counts/counts.tsv`, `manifest/manifest.tsv`,
`markers/specificity_report.tsv`, `markers/mask_intervals.tsv`,
`emit/all_cds.ffn` (**35 GB** — cross-map target universe), `emit/reps.ffn`
(189 MB), `low_marker/skani.sparse.tsv` (3.1M rows).

---

## 4. Roadmap

### Phase 0 — Diagnostics ✅
Leakage extraction, radial visualisation, lower-bound phylo merge probe, and the
identity/locality characterisation above. Done this session.

### Phase 1 — Medoid selection ⬜ (local)
- `bin/clade_medoid.py`: given a clade (species list) + `skani.sparse.tsv` +
  `manifest.tsv`, return the medoid genome (max mean ANI; handle missing pairs in
  the sparse table by treating them as below-screen / low ANI).
- Emit an audit row per clade: medoid accession, mean/min ANI, n genomes,
  n missing pairs (sparsity warning).

### Phase 2 — Merge re-evaluation process ⬜ (M3 / Nextflow)
- New subworkflow `modules/local/merge_reeval.nf`, cloning the relax path:
  - `MERGE_EMIT`: prevalence candidates for the merged clade (reuse
    `phylo_merge`/`merge_gain` candidate logic) → medoid CDS → query FASTA +
    clade relabel map.
  - `MERGE_CROSSMAP`: `mmseqs easy-search` of medoid CDS vs `all_cds.ffn`
    **excluding in-clade genome idxs** (build an exclusion set from the manifest).
  - `MERGE_GUARD`: `specificity_guard.py` with masking → clean marker set +
    per-marker mask report for the merged clade.
- Output per clade: `merged_markers`, delta vs baseline, recovered-by-masking
  count, and the fresh cross-map footprint.

### Phase 2b — Core-gene mask rescue ✅ IMPLEMENTED
Orthogonal to merging: recover markers from the *locally divergent window* of an
otherwise-shared core gene. For flagged species, drop the uniqueness gate (b),
cross-map **every** core gene against all-CDS, and keep the ones the guard
recovers by masking. Rationale: overall uniqueness (`out_prev`) and local
maskability are orthogonal, so gate (b) filters on the wrong axis and discards
genes whose divergent window is species-diagnostic.

- Selector [`bin/select_core.py`](../bin/select_core.py): all core genes
  (`in_prev ≥ min_in`) for flagged species, uniqueness dropped, already-emitted
  excluded. Verified on the v2 run: **2,780** candidates for ovatus (= 2862 core −
  82 emitted), 62,769 across 23 flagged species.
- Subworkflow [`modules/local/core_rescue.nf`](../modules/local/core_rescue.nf):
  `SELECT_CORE → CORE_EMIT → CORE_CROSSMAP → CORE_GUARD → SELECT_CORE_RESCUED →
  MERGE_CORE`. `CORE_GUARD` runs `specificity_guard.py` **with** `--target
  all_cds` (mask recovery on), unlike the clean-only relax guard.
- Wired into [`workflows/markers.nf`](../workflows/markers.nf) (step 9d), chains
  onto `markers_final` after the relax path. Params: `core_rescue` (default
  `false`, opt-in — cross-mapping ~63k core genes is mmseqs-heavy),
  `core_rescue_max_per_clade` (0 = every core gene). `nextflow lint` clean;
  `-preview` DAG resolves.
- **Within-species window conservation** (roadmap §5 concern) is covered by the
  `in_prev ≥ 0.80` core gate + keeping the whole gene, not just the window; the
  finer "is the *clean window itself* conserved across the species' genomes" check
  is a future refinement.
- Not yet run (needs mmseqs + the 35 GB all_cds on M3). To run:
  `nextflow run ... --core_rescue true`.

### Phase 3 — Iterative driver ⬜
- Walk the phylogeny per flagged species (as `phylo_merge` does), calling Phase 2
  at each node, until a **stop criterion** is met (see §5).
- Record the trajectory: clade size, leakage% removed, markers recovered (real,
  not lower bound), diversity cost (n species collapsed).
- Choose the **smallest clade** that clears the marker bar with the **cleanest ANI
  support** (see gate in §5).

### Phase 4 — Merge decision + DB integration ⬜
- Apply the ANI gate: only commit a merge where the clade has **no clean within/
  between ANI gap** (grey-zone) — never collapse cleanly-separable taxa.
- Relabel committed merges in the marker DB as a named complex
  (e.g. `Bacteroides ovatus/xylanisolvens complex`) and report classification at
  that resolution; keep an explicit provenance map (merged GTDB species → complex).
- Regenerate the low-marker report to show pre/post-merge marker counts.

### Phase 5 — Allelic-marker fallback ⬜ (out of scope for merging, tracked here)
For complexes where **no** merge recovers markers (leakage 100 % identical across
the whole complex), gene-presence markers cannot work. Track a separate track:
allelic/SNP markers on the shared core (StrainPhlAn/metaSNV-style) — the stratum
`out_prev` currently discards. Decision needed on whether species-level resolution
of grey-zone complexes is a project goal at all.

---

## 5. Decision points / open questions

- **Representative choice.** Medoid genome (proposed) vs union of member reps.
  Medoid is cheaper and central; union is more complete but re-introduces
  redundancy. Default: medoid + per-cluster fallback.
- **Stop criterion for the walk.** Options: (a) first clade reaching
  `low_marker_threshold` markers; (b) first clade where marker gain plateaus;
  (c) leakage-removed ≥ X %. Recommend (a), tie-broken by smallest clade / best ANI.
- **ANI gate.** A merge must be *authorised* by genome-level similarity, not by
  marker deficit. Only merge clades with overlapping within/between ANI (no clean
  gap). This keeps merges monophyletic **and** evidence-based — never merge to
  chase a count.
- **Diversity cost budget.** Because ovatus's leakage is polyphyletic, the clade
  that internalises it may be large. Define a max acceptable number of species
  collapsed per recovered marker, below which we prefer Phase 5 (allelic) or
  accept low resolution.
- **Reporting.** How merged taxa surface in `marker_db.tsv` and downstream
  classification; provenance so a merge is reversible/auditable.

---

## 6. Risks

- **May still recover little for ovatus.** If the shared-core genes are also
  ~100 % identical to out-of-clade species, even a correct re-emit drops them.
  The re-evaluation tells us the truth; it does not guarantee recovery. Value is
  the same for species with **monophyletic** leakage, which should recover well.
- **skani sparsity.** The medoid may be computed from an incomplete ANI matrix;
  flag low-coverage clades.
- **Cost.** Fresh cross-map per candidate clade per flagged species is
  `mmseqs`-heavy; cache by clade signature and only walk flagged species.

---

## 7. Success criteria

1. For each flagged species, a **true** (not lower-bound) merged-marker count at
   each phylogenetic node, with the diversity cost.
2. A clear split of flagged species into: **recoverable by monophyletic merge**
   (markers return, ANI gate passes) vs **allelic-only** (no merge recovers).
3. `B. ovatus` resolved to a defensible outcome: either a named complex with a
   real marker set, or an explicit "not gene-content-separable" verdict backed by
   the re-emit evidence.
