# meta-mage — memory-scaling ROADMAP (8k → 800k genomes)

Status: **draft / planning.** This document audits every stage for memory behaviour
at 800k genomes and lays out the internal changes needed to get there. It is about
*memory efficiency*, not features. Nothing here is done yet.

## The target

| Quantity | 8k run (now) | 800k run (goal) | Factor |
|---|---|---|---|
| Genomes | ~8k | ~800k (≈ all GTDB r220) | 100× |
| Proteins (`.faa`) | ~28M | ~2.8B | 100× |
| Species reps | ~hundreds | ~110k | — |
| Rep CDS (`reps.ffn`) | ~1M | ~385M (~0.4–1 TB) | — |
| Cluster TSV rows | ~28M | ~2.8B (~150 GB ×2) | 100× |
| Markers (cap 200/clade) | ~thousands | up to ~30M, many copies each | — |
| Scatter tasks (PYRODIGAL+REHEADER×2) | ~24k | ~2.4M | 100× |

Cluster ceiling (`m3c`): **1843 GB RAM, 128 cores, 14-day walltime**, `cpu2-hm`
high-mem queue. Generous — but three stages below blow past *even that*.

## The invariant we want

> **Hold nothing in RAM that scales with the protein/CDS corpus.** Only
> `idx → lineage` (≈ #genomes, compact) and per-clade counters (≈ taxonomy size)
> may be resident. Everything else must **stream a sorted/grouped input one block
> at a time**, and sequence payloads must be **random-accessed via an index**
> (faidx), never fully loaded.

`bin/aggregate_counts.py` already embodies this (streams the rep-sorted cluster
TSV, holds only the current cluster's genome set + the interned `idx→lineage`
map). The roadmap is largely: **make the other stages behave like COUNTS.**

## Stage-by-stage audit

Peak RAM is the *current-design* estimate at 800k. ✅ = already scales.

| Stage | File | Holds in RAM now | ~Peak @800k | Verdict |
|---|---|---|---|---|
| BUILD_MANIFEST | `bin/build_manifest.py` | whole sheet as a `rows` list | 1–2 GB | ⚠️ wasteful, survivable |
| PYRODIGAL / REHEADER | modules | per-task only | small | ⚠️ **2.4M tasks** (filesystem/scheduler wall) |
| MMSEQS_CREATEDB | nf-core | — | low | ✅ |
| MMSEQS_LINCLUST ×2 | nf-core | k-mer index | **unbounded w/o cap** | ⚠️ needs `--split-memory-limit` |
| MMSEQS_CREATETSV | nf-core | — | low RAM / 150 GB disk | ✅ (disk) |
| COUNTS | `bin/aggregate_counts.py` | `idx→lineage` + current cluster set | ~0.3–0.5 GB | ✅ |
| SCORE | `bin/score_markers.py` | **all passing candidates** (`by_clade`) | 1–5 GB, unbounded in principle | ⚠️ bound it |
| EMIT_REPS | `bin/emit_reps.py` | **`load_fasta(reps.ffn)` — every rep CDS** | **0.4–1 TB** | ❌ **OOM** |
| CROSSMAP | `modules/local/crossmap.nf` | mmseqs easy-search over 2.8B CDS target | **unbounded w/o cap** | ♻️ **redesign** → neighborhood read-sim, sharded (see below) |
| SPECIFICITY | `bin/specificity_guard.py` | `idx→lineage` + `offenders` (per marker) | 0.3–2 GB | ⚠️ monitor |
| DIVERSITY | `bin/diversity.py` | **`read_groups` — the entire marker FASTA** | **100s of GB** | ❌ **OOM** |
| EMIT_DB | `bin/emit_db.py` | **`seen` — every emitted sequence** | **100s of GB** | ❌ **OOM** |

### The three OOM blockers (scale with the sequence corpus, not the taxonomy)

- **EMIT_REPS** — `emit_reps.py:101` `rep_seqs = load_fasta(args.reps_ffn)` slurps
  all ~385M rep CDS into a `dict[str,str]`. The code already flags this
  (`emit_reps.py:99-100`). Only a tiny, capped subset of sequences is ever
  written, so a full load is almost pure waste.
- **DIVERSITY** — `diversity.py:80` `read_groups(args.fasta)` loads the whole
  marker FASTA into `groups`, then builds k-mer sets for all of it.
- **EMIT_DB** — `emit_db.py:60,71` `seen[key].add(seq)` keeps every emitted
  sequence string resident to dedup exact copies.

All three have the same cure: **stream a marker-grouped FASTA and hold only the
current marker's copies**, plus (for EMIT_REPS) fetch payloads by index.

## Cross-cutting enablers (do these first — they unlock the stage fixes)

1. **Canonical marker ordering `rank|clade|cluster`.** If `EMIT_REPS` writes its
   FASTA already grouped by `(rank, clade, cluster)`, then DIVERSITY and EMIT_DB
   can each process one marker block at a time and drop their global maps. One
   ordering contract fixes two OOMs.
2. **faidx on `reps.ffn`.** Build a `.fai` (samtools faidx, or a stdlib
   byte-offset index) so EMIT_REPS random-accesses only the member sequences it
   emits. Turns a 0.4–1 TB load into O(payload written).
3. **Sorted `counts.tsv` by `(rank, clade)`.** Lets SCORE keep a bounded top-k
   heap per clade instead of `by_clade` holding all passing candidates. The sort
   is external/on-disk (cheap; COUNTS already relies on `LC_ALL=C sort`).
4. **mmseqs `--split-memory-limit`.** Add to LINCLUST `ext.args`, sized to the
   allocation (e.g. `--split-memory-limit $((task.memory*0.8))`). Bounds RAM and
   spills to disk; this is how mmseqs clusters billions of seqs. (CROSSMAP no
   longer needs this — it is being replaced; see the redesign below.)

## Workstreams (prioritised)

### P0 — OOM blockers (pipeline cannot finish 800k without these)
- [ ] **EMIT_REPS**: faidx `reps.ffn` + fetch-on-demand; emit grouped by
      `rank|clade|cluster`; replace `idx_lineage` dict-of-dicts with a tuple
      (as COUNTS/SPECIFICITY already do).
- [ ] **DIVERSITY**: stream the grouped FASTA one marker block at a time; hold
      only that block's copies. Optionally switch to `mash sketch -i`/`mash dist`
      per block (the module now has a mash container path if we go that way).
- [ ] **EMIT_DB**: stream the grouped FASTA; dedup within the current marker
      block only (bounded `seen`).

### P1 — scale/throughput walls
- [ ] **Chunk the per-genome scatter.** Batch G genomes per PYRODIGAL/REHEADER
      task (channel `.collate(G)` or a manifest-shard input) to cut ~2.4M tasks
      to ~24k. Job arrays (already added) help submission; chunking cuts the
      work-dir/staging explosion that arrays don't.
- [ ] **mmseqs `--split-memory-limit`** on LINCLUST; ensure large, fast scratch
      (`-work-dir` on Lustre is fine; confirm `TMPDIR`).
- [ ] **CROSSMAP/specificity redesign** — replace the global all-CDS mmseqs
      search with a neighborhood-scoped read-simulation guard, sharded by
      neighborhood. Full design in the section below. Removes the single heaviest
      step and its memory concern entirely.

### P2 — bounded but worth tightening
- [ ] **SCORE**: sorted-counts streaming + per-clade top-k heap (needs enabler 3).
- [ ] **BUILD_MANIFEST**: two-pass (count species, then stream-emit) so the full
      `rows` list is never resident. Low urgency (~1–2 GB).
- [ ] **SPECIFICITY**: monitor `offenders`; if it grows with markers, key the
      report by streaming the marker table in the canonical order and joining
      against a sorted hits file rather than a dict.

## CROSSMAP / specificity — redesign (decided 2026-07-14)

Replace the global mmseqs `easy-search` guard (`crossmap.nf` + `specificity_guard.py`)
with a **neighborhood-scoped, read-simulation** test that measures the actual
classifier failure mode directly: *reads originating outside a marker's clade that
accumulate on that marker.*

**Rationale.** The global out-of-clade filter already runs in protein space
(SCORE `out_prevalence` over all genomes), so the nucleotide residual is *local*
(sister taxa within genus/family). And what harms the classifier is *read
mis-assignment*, not pident per se — so simulate reads and count where they land,
instead of thresholding alignment identity.

**Algorithm — per neighborhood `P`** (P = the marker's ancestor at
`rank + crossmap_neighborhood_depth`; e.g. depth 2 → species markers judged within
their family):
1. **Simulate/tile reads** (length = the downstream classifier's read length) from
   every genome in `P`.
2. **Map** them to the marker set of `P`.
3. Per marker, `off_frac = reads_from_outside_marker_clade / reads_on_marker`
   (clade evaluated at the marker's own rank; a genome's reads mapping to a marker
   of a *different* member of `P` is the bad signal).
4. **Flag/drop** if `off_frac > max_offspecies_read_frac`, gated on a minimum read
   count so ultra-low-coverage markers aren't judged on noise.

**Why it scales.** Naturally sharded by `P`: target = only `P`'s genomes + markers.
Work ≈ `Σ_P (reads_P × markers_P)` — bounded by neighborhood size, not by N. Each
shard fits a `cpu1` node and shards run in parallel. This *replaces*
`--split-memory-limit` for this stage.

**New params.**
- `crossmap_neighborhood_depth = 2` — ancestor levels up defining the off-target
  universe (1 = parent/genus, 2 = grandparent/family). Configurable per run.
- `crossmap_read_len` — simulated read length. **Must match the downstream
  classifier's**, or the guard tests the wrong mapping regime.
- `max_offspecies_read_frac` — per-marker flag threshold.
- `crossmap_min_reads` — minimum reads on a marker before it is judged.

**Tooling.** Deterministic tiling of CDS into `crossmap_read_len` windows is
stdlib-Python (reproducible, no RNG) and is the conservative floor; a real
simulator (`art`/`dwgsim`) adds sequencing-error realism (*more* sensitive to
cross-mapping) if wanted. Mapping needs one short-read-mapper container — ideally
the **same mapper + settings the classifier uses**; else minimap2 `-sr` / bowtie2.
(Same pattern as the mash container just added to DIVERSITY.)

**Open coupling to resolve before implementing.** The guard is only meaningful if
its read length and mapper mirror the tool that will consume this DB. Confirm the
downstream classifier's mapping regime and pin `crossmap_read_len` + mapper to it.

**Residual risk (accepted).** Cross-*neighborhood* near-identity (e.g. recent HGT
shared at high nt identity across families that protein out-prevalence missed).
Depth 2 (family) is the safety margin; a mobile-element screen could mop up the
rest if it ever proves necessary.

## Validation strategy

- Instrument every `bin/*.py` process with `/usr/bin/time -v` (peak RSS) via a
  `beforeScript`/wrapper; capture in the trace. Peak RSS per stage is the metric.
- Build a **synthetic scaling harness**: replicate the Bacteroides sheet ×N
  (10k, 50k, 200k) to plot peak RSS vs. genomes per stage *before* committing a
  full 800k run. A stage that scales must stay flat or grow with taxonomy, not
  with N.
- Track intermediate sizes (`all_proteins.faa`, cluster TSVs, marker FASTA) so
  disk is provisioned, not discovered at 90%.
- Gate: no stage may exceed a set fraction (e.g. 25%) of a `cpu1` node before we
  reach for `cpu2-hm`.

## Not memory, but will bite at 800k (track separately)
- Disk: `all_proteins.faa` ~1 TB, cluster TSVs ~150 GB ×2, marker FASTA large.
  Compress intermediates; prefer `collectFile(..., sort: false)` where ordering
  is imposed later.
- Walltime: LINCLUST ×2 and CROSSMAP over billions of seqs — size against the
  14-day `cpu3-long` queue; rely on the retry-with-more (`base.config`) ladder.
- Nextflow overhead: ~2.4M tasks strains the executor and `.nextflow` history;
  chunking (P1) is the mitigation.
