#!/usr/bin/env python3
"""
Build a self-contained HTML report of the marker-discovery run.

The report reconstructs, per species, the funnel that turns a pangenome into a
set of clade-specific markers, and records at every stage how many genes were
removed and why:

    pangenome  ->  core  ->  clade-specific candidate  ->  selected  ->  final
       |            |               |                         |           |
   all gene     present in     specific enough           within the   survived the
   families     >= core_prev    (out-prevalence <=        per-clade    nucleotide
   in the       of the         max_out, in-prevalence     marker cap   cross-map
   species      species        >= min_in, clade big       (max_per_    guard (QC)
                genomes         enough)                    clade)

Everything is derived from the pipeline's own tables so the numbers match what
the DB actually contains:
  * manifest.tsv        -- kept genomes, lineage, per-species genome counts (N)
  * dropped_species.tsv -- species removed by the genome-count prefilter
  * counts.tsv          -- (cluster, rank, clade, in_count, marker_total)
  * clade_sizes.tsv     -- (rank, clade, size) + __TOTAL__
  * markers.tsv         -- SCORE output: the selected (post-cap) candidates
  * specificity_report  -- per-marker cross-map QC verdict (optional)

Species-rank clusters are tagged "S:" in counts.tsv; that is the vocabulary the
species funnel is built from. The QC-removal and per-genus sections cover every
rank. Output is a single dependency-free HTML file (stdlib only).
"""
import argparse
import html
import json
import os
from collections import defaultdict
from datetime import datetime, timezone

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


# --------------------------------------------------------------------------- IO
def load_manifest(path):
    """Return (species_size, species_genus, n_total).

    species_size[species]  = # genomes of that species (kept after prefilter)
    species_genus[species] = the genus the species belongs to
    n_total                = total kept genomes (out-clade denominator base)
    """
    species_size = defaultdict(int)
    species_genus = {}
    n_total = 0
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sp = f[col["species"]]
            n_total += 1
            if sp == "NA":
                continue
            species_size[sp] += 1
            species_genus[sp] = f[col["genus"]]
    return species_size, species_genus, n_total


def load_dropped(path):
    """species -> (genus, n_genomes, min_required) for prefilter-dropped species."""
    dropped = {}
    if not path or not os.path.exists(path):
        return dropped
    with open(path) as fh:
        header = fh.readline()
        if not header:
            return dropped
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            dropped[f[0]] = (f[1], int(f[2]), int(f[3]))
    return dropped


def load_specificity(path):
    """(rank, clade, cluster) -> dict verdict. Empty when the guard didn't run."""
    verdicts = {}
    if not path or os.path.basename(path) == "NO_FILE" or not os.path.exists(path):
        return verdicts
    with open(path) as fh:
        header = fh.readline()
        if not header:
            return verdicts
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            rank, clade, cluster, passed, n_off, max_id, worst = f[:7]
            verdicts[(rank, clade, cluster)] = {
                "pass": passed == "1",
                "n_offtarget": int(n_off),
                "max_id": float(max_id),
                "worst": worst,
            }
    return verdicts


def load_selected(path):
    """(rank, clade, cluster) -> score for SCORE's selected (post-cap) markers."""
    selected = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            selected[(f[0], f[1], f[2])] = f[7] if len(f) > 7 else ""
    return selected


# ------------------------------------------------------------------- analysis
def build_species_funnel(counts_path, species_size, n_total, p):
    """Stream counts.tsv (species rows) into a per-species funnel.

    Returns species -> counters:
        pangenome, core, non_core, specific, non_specific, below_min_in,
        clade_too_small (0/pangenome), plus the set of specific candidate keys.
    """
    funnel = defaultdict(lambda: {
        "pangenome": 0, "core": 0, "non_core": 0,
        "specific": 0, "non_specific": 0, "below_min_in": 0,
        "clade_too_small": False, "specific_keys": [],
    })
    with open(counts_path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            cluster, rank, clade, in_count, marker_total = f[:5]
            if rank != "species" or clade == "NA":
                continue
            size = species_size.get(clade)
            if not size:
                continue
            in_count = int(in_count)
            marker_total = int(marker_total)
            fn = funnel[clade]
            fn["pangenome"] += 1
            in_prev = in_count / size
            if in_prev < p["core_prevalence"]:
                fn["non_core"] += 1
                continue
            fn["core"] += 1
            # Recompute SCORE's specificity filter for this core gene.
            if size < p["min_clade_size"]:
                fn["clade_too_small"] = True
                fn["non_specific"] += 1
                continue
            out_denom = n_total - size
            out_prev = (marker_total - in_count) / out_denom if out_denom > 0 else 1.0
            if in_prev < p["min_in"]:
                fn["below_min_in"] += 1
            elif out_prev <= p["max_out"]:
                fn["specific"] += 1
                fn["specific_keys"].append((rank, clade, cluster))
            else:
                fn["non_specific"] += 1
    return funnel


# ---------------------------------------------------------------------- render
def esc(x):
    return html.escape(str(x))


def pct(n, d):
    return f"{100.0 * n / d:.1f}%" if d else "-"


def card(label, value, sub=""):
    sub = f"<div class='sub'>{esc(sub)}</div>" if sub else ""
    return (f"<div class='card'><div class='val'>{esc(value)}</div>"
            f"<div class='lbl'>{esc(label)}</div>{sub}</div>")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--clade_sizes", required=True)  # accepted for provenance
    ap.add_argument("--markers", required=True, help="SCORE selected markers")
    ap.add_argument("--dropped", help="dropped_species.tsv from the prefilter")
    ap.add_argument("--specificity", help="specificity_report.tsv (or NO_FILE)")
    ap.add_argument("--min_in", type=float, required=True)
    ap.add_argument("--max_out", type=float, required=True)
    ap.add_argument("--min_clade_size", type=int, required=True)
    ap.add_argument("--max_per_clade", type=int, required=True)
    ap.add_argument("--core_prevalence", type=float, default=0.90)
    ap.add_argument("--min_genomes_per_species", type=int, default=1)
    ap.add_argument("--max_table_rows", type=int, default=5000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    p = {
        "min_in": args.min_in, "max_out": args.max_out,
        "min_clade_size": args.min_clade_size, "max_per_clade": args.max_per_clade,
        "core_prevalence": args.core_prevalence,
    }

    species_size, species_genus, n_total = load_manifest(args.manifest)
    dropped = load_dropped(args.dropped)
    verdicts = load_specificity(args.specificity)
    selected = load_selected(args.markers)
    qc_ran = len(verdicts) > 0
    funnel = build_species_funnel(args.counts, species_size, n_total, p)

    # ---- per-species selected / capped / final, from the actual marker tables.
    sp_selected = defaultdict(int)   # species -> # selected (post-cap) markers
    sp_final = defaultdict(int)      # species -> # markers surviving QC
    for (rank, clade, cluster) in selected:
        if rank != "species":
            continue
        sp_selected[clade] += 1
        v = verdicts.get((rank, clade, cluster))
        if not qc_ran or (v and v["pass"]):
            sp_final[clade] += 1

    # ---- markers-by-rank summary (all ranks).
    rank_selected = defaultdict(int)
    rank_final = defaultdict(int)
    rank_qc_dropped = defaultdict(int)
    for (rank, clade, cluster) in selected:
        rank_selected[rank] += 1
        v = verdicts.get((rank, clade, cluster))
        if qc_ran and v and not v["pass"]:
            rank_qc_dropped[rank] += 1
        else:
            rank_final[rank] += 1

    # ---- QC removals (cross-map guard) detail, all ranks.
    qc_removed = [
        (rank, clade, cluster, v["n_offtarget"], v["max_id"], v["worst"])
        for (rank, clade, cluster), v in verdicts.items() if not v["pass"]
    ]
    qc_removed.sort(key=lambda r: (-r[3], -r[4]))

    # ---- per-genus removal breakdown.
    # total species per genus = kept (in manifest) + prefilter-dropped.
    genus_species = defaultdict(set)
    for sp, g in species_genus.items():
        genus_species[g].add(sp)
    genus_prefilter = defaultdict(set)
    for sp, (g, _n, _m) in dropped.items():
        genus_prefilter[g].add(sp)
        genus_species[g].add(sp)  # dropped species aren't in the manifest
    zero_marker_species = sorted(
        sp for sp in species_size
        if sp_final.get(sp, 0) == 0
    )
    genus_zero = defaultdict(set)
    for sp in zero_marker_species:
        genus_zero[species_genus.get(sp, "NA")].add(sp)

    # ---- totals.
    tot_pangenome = sum(fn["pangenome"] for fn in funnel.values())
    tot_core = sum(fn["core"] for fn in funnel.values())
    tot_specific = sum(fn["specific"] for fn in funnel.values())
    tot_selected = sum(sp_selected.values())
    tot_final = sum(sp_final.values())
    n_species = len(species_size)
    n_genera = len({g for g in species_genus.values() if g != "NA"})

    html_out = render_html(
        args, p, qc_ran, n_total, n_species, n_genera,
        species_size, species_genus, funnel, sp_selected, sp_final,
        rank_selected, rank_final, rank_qc_dropped,
        qc_removed, dropped, genus_species, genus_prefilter, genus_zero,
        zero_marker_species,
        totals=dict(pangenome=tot_pangenome, core=tot_core, specific=tot_specific,
                    selected=tot_selected, final=tot_final),
    )
    with open(args.out, "w") as fh:
        fh.write(html_out)
    print(f"build_report: wrote {args.out} "
          f"({n_species} species, {tot_final} final markers, "
          f"{len(qc_removed)} removed by QC)")


def render_html(args, p, qc_ran, n_total, n_species, n_genera,
                species_size, species_genus, funnel, sp_selected, sp_final,
                rank_selected, rank_final, rank_qc_dropped,
                qc_removed, dropped, genus_species, genus_prefilter, genus_zero,
                zero_marker_species, totals):
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # Summary cards.
    cards = "".join([
        card("Genomes", f"{n_total:,}"),
        card("Species", f"{n_species:,}"),
        card("Genera", f"{n_genera:,}"),
        card("Gene families", f"{totals['pangenome']:,}", "species-rank, summed"),
        card("Core genes", f"{totals['core']:,}", f">= {p['core_prevalence']:.0%} prevalence"),
        card("Specific candidates", f"{totals['specific']:,}", "pre-cap"),
        card("Selected markers", f"{totals['selected']:,}", f"cap {p['max_per_clade']}/clade"),
        card("Final markers", f"{totals['final']:,}",
             "after cross-map QC" if qc_ran else "QC not run"),
    ])

    # Species funnel table.
    species_rows = sorted(
        set(funnel) | set(species_size),
        key=lambda s: (species_genus.get(s, "NA"), s),
    )
    truncated = len(species_rows) > args.max_table_rows
    body_rows = []
    for sp in species_rows[: args.max_table_rows]:
        fn = funnel.get(sp, {})
        pan = fn.get("pangenome", 0)
        core = fn.get("core", 0)
        spec = fn.get("specific", 0)
        sel = sp_selected.get(sp, 0)
        fin = sp_final.get(sp, 0)
        capped = max(0, spec - sel)
        qc_drop = max(0, sel - fin) if qc_ran else 0
        cls = " class='warn'" if fin == 0 else ""
        body_rows.append(
            f"<tr{cls}><td class='l'>{esc(sp)}</td>"
            f"<td class='l dim'>{esc(species_genus.get(sp, 'NA'))}</td>"
            f"<td>{species_size.get(sp, 0)}</td>"
            f"<td>{pan}</td><td>{core}</td>"
            f"<td>{fn.get('non_core', 0)}</td>"
            f"<td>{spec}</td>"
            f"<td>{fn.get('non_specific', 0) + fn.get('below_min_in', 0)}</td>"
            f"<td>{capped}</td><td>{qc_drop}</td>"
            f"<td class='strong'>{fin}</td></tr>"
        )
    trunc_note = (f"<p class='note'>Showing the first {args.max_table_rows:,} of "
                  f"{len(species_rows):,} species (sorted by genus). Full per-marker "
                  f"data is in the pipeline's counts/markers TSVs.</p>"
                  if truncated else "")

    # Markers-by-rank table.
    rank_body = []
    for rank in RANKS:
        if rank_selected.get(rank, 0) == 0:
            continue
        rank_body.append(
            f"<tr><td class='l'>{esc(rank)}</td>"
            f"<td>{rank_selected.get(rank, 0)}</td>"
            f"<td>{rank_qc_dropped.get(rank, 0) if qc_ran else '-'}</td>"
            f"<td class='strong'>{rank_final.get(rank, 0)}</td></tr>"
        )

    # QC-removal detail.
    if not qc_ran:
        qc_section = ("<p class='note'>The nucleotide cross-map guard was not run "
                      "(<code>--specificity false</code>), so no markers were removed "
                      "by QC.</p>")
    elif not qc_removed:
        qc_section = ("<p class='note ok'>No markers were removed by the cross-map "
                      "guard — every selected marker is nucleotide-unique to its "
                      "clade.</p>")
    else:
        shown = qc_removed[: args.max_table_rows]
        qc_body = "".join(
            f"<tr><td class='l dim'>{esc(rank)}</td><td class='l'>{esc(clade)}</td>"
            f"<td class='l dim'>{esc(cluster)}</td><td>{n_off}</td>"
            f"<td>{max_id:.1%}</td><td class='l'>{esc(worst)}</td></tr>"
            for (rank, clade, cluster, n_off, max_id, worst) in shown
        )
        more = (f"<p class='note'>Showing {len(shown):,} of {len(qc_removed):,} "
                f"removed markers.</p>" if len(qc_removed) > len(shown) else "")
        qc_section = (
            f"<p class='note'>These selected markers were dropped because their "
            f"nucleotide sequence also occurs in a genome <em>outside</em> the target "
            f"clade at &ge; {args.min_in:.0%}-scale identity — such a marker would "
            f"steal reads from the off-target taxon. Removal is triggered per "
            f"(rank, clade, cluster).</p>"
            "<table><thead><tr><th>Rank</th><th>Clade</th><th>Cluster</th>"
            "<th>Off-target hits</th><th>Max identity</th>"
            "<th>Worst off-target clade</th></tr></thead>"
            f"<tbody>{qc_body}</tbody></table>{more}"
        )

    # Removed-species: per-genus breakdown.
    genus_rows = []
    genera_sorted = sorted(
        genus_species,
        key=lambda g: (-len(genus_prefilter.get(g, ())), g),
    )
    for g in genera_sorted:
        if g == "NA":
            continue
        total_sp = len(genus_species[g])
        n_pre = len(genus_prefilter.get(g, ()))
        n_zero = len(genus_zero.get(g, ()))
        if n_pre == 0 and n_zero == 0:
            continue
        genus_rows.append(
            f"<tr><td class='l'>{esc(g)}</td><td>{total_sp}</td>"
            f"<td>{n_pre}</td><td>{pct(n_pre, total_sp)}</td>"
            f"<td>{n_zero}</td><td>{pct(n_zero, total_sp)}</td></tr>"
        )
    if genus_rows:
        genus_table = (
            "<table><thead><tr><th>Genus</th><th>Species (total)</th>"
            "<th>Removed: too few genomes</th><th>%</th>"
            "<th>Kept but 0 markers</th><th>%</th></tr></thead>"
            f"<tbody>{''.join(genus_rows)}</tbody></table>"
        )
    else:
        genus_table = "<p class='note ok'>No species were removed.</p>"

    # Prefilter drop list.
    if dropped:
        pre_body = "".join(
            f"<tr><td class='l'>{esc(sp)}</td><td class='l dim'>{esc(g)}</td>"
            f"<td>{n}</td><td>{req}</td></tr>"
            for sp, (g, n, req) in sorted(dropped.items())
        )
        pre_section = (
            f"<p class='note'>Dropped before indexing by "
            f"<code>--min_genomes_per_species={args.min_genomes_per_species}</code>: "
            f"a species with too few genomes gives meaningless within-species "
            f"prevalence and pollutes higher-rank clade sizes.</p>"
            "<table><thead><tr><th>Species</th><th>Genus</th>"
            "<th>Genomes</th><th>Required</th></tr></thead>"
            f"<tbody>{pre_body}</tbody></table>"
        )
    else:
        pre_section = ("<p class='note ok'>No species were removed by the "
                       "genome-count prefilter.</p>")

    # Params block.
    param_items = "".join(
        f"<span class='pill'>{esc(k)} = {esc(v)}</span>" for k, v in [
            ("core_prevalence", p["core_prevalence"]),
            ("min_in_prevalence", args.min_in),
            ("max_out_prevalence", args.max_out),
            ("min_clade_size", args.min_clade_size),
            ("max_markers_per_clade", args.max_per_clade),
            ("min_genomes_per_species", args.min_genomes_per_species),
            ("cross_map_QC", "on" if qc_ran else "off"),
        ]
    )

    # Column-header tooltips: definition + exact filter + the run's parameters.
    cp = f"{p['core_prevalence']:.0%}"
    mi, mo = args.min_in, args.max_out
    mcs, mpc = args.min_clade_size, args.max_per_clade
    mgps = args.min_genomes_per_species
    qid = f"{args.min_in:.0%}"
    tips = {
        "genomes": (
            "<b>Genomes</b>"
            "<div>Genomes of this species kept after the genome-count prefilter. "
            "This is the denominator N used for every in-species prevalence.</div>"
            f"<div class='p'>min_genomes_per_species = {mgps}</div>"
        ),
        "pangenome": (
            "<b>Pangenome</b>"
            "<div>Every gene family (cluster) seen in &ge;1 genome of the "
            "species &mdash; the starting pool of the funnel.</div>"
            "<div class='f'>Filter: none</div>"
            "<div class='p'>= count of all species-rank clusters for this "
            "species in counts.tsv</div>"
        ),
        "core": (
            "<b>Core</b>"
            "<div>Gene families present in most of the species' genomes.</div>"
            f"<div class='f'>keep if&nbsp;&nbsp;in_prev &ge; {cp}</div>"
            "<div class='p'>in_prev = in_count / genomes"
            f"<br>parameter: core_prevalence = {p['core_prevalence']}</div>"
        ),
        "not_core": (
            "<b>Not core</b>"
            "<div>Gene families dropped at the core step.</div>"
            f"<div class='f'>in_prev &lt; {cp}</div>"
            "<div class='p'>= Pangenome &minus; Core</div>"
        ),
        "specific": (
            "<b>Specific</b>"
            "<div>Core genes clade-specific enough to be markers. "
            "All three must hold:</div>"
            f"<div class='f'>genomes &ge; {mcs}<br>"
            f"in_prev &ge; {mi}<br>"
            f"out_prev &le; {mo}</div>"
            "<div class='p'>out_prev = (marker_total &minus; in_count) / "
            "(total_genomes &minus; genomes)"
            f"<br>parameters: min_clade_size = {mcs}, min_in = {mi}, "
            f"max_out = {mo}</div>"
        ),
        "not_specific": (
            "<b>Not specific</b>"
            "<div>Core genes that failed the specificity test for any reason:</div>"
            f"<div class='f'>genomes &lt; {mcs}"
            f"<br>&nbsp;&nbsp;OR&nbsp; in_prev &lt; {mi}"
            f"<br>&nbsp;&nbsp;OR&nbsp; out_prev &gt; {mo}</div>"
            "<div class='p'>= Core &minus; Specific</div>"
        ),
        "capped": (
            "<b>Capped</b>"
            "<div>Specific candidates discarded because the per-clade marker cap "
            "was reached; SCORE keeps only the top-scoring markers per clade.</div>"
            f"<div class='f'>rank &gt; {mpc} by score</div>"
            f"<div class='p'>parameter: max_per_clade = {mpc}"
            "<br>= Specific &minus; Selected</div>"
        ),
        "selected": (
            "<b>Selected</b>"
            "<div>Specific candidates kept after the per-clade cap &mdash; the "
            "markers SCORE actually selected (pre-QC).</div>"
            f"<div class='f'>top &le; {mpc} per clade, by score</div>"
            f"<div class='p'>parameter: max_per_clade = {mpc}</div>"
        ),
        "qc_dropped": (
            "<b>QC dropped</b>"
            "<div>Selected markers removed by the nucleotide cross-map guard: the "
            "marker sequence also maps to a genome <em>outside</em> the target "
            f"clade at &ge; {qid}-scale identity, so it would steal reads from the "
            "off-target taxon.</div>"
            "<div class='f'>specificity_report verdict = fail</div>"
            "<div class='p'>= Selected &minus; Final"
            "<br><b>not</b> included in Final</div>"
        ),
        "final": (
            "<b>Final</b>"
            "<div>Markers that survived every stage, including cross-map QC. "
            "This is the marker set written to the DB.</div>"
            "<div class='f'>= Selected &minus; QC dropped</div>"
            "<div class='p'>QC-dropped genes are excluded from this count.</div>"
        ),
    }
    tips_json = json.dumps(tips)

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>meta-mage marker report</title>
<style>
:root {{
  --bg:#0f1115; --panel:#171a21; --panel2:#1d212b; --fg:#e6e9ef; --dim:#9aa4b2;
  --line:#2a2f3a; --accent:#5aa9e6; --warn:#e0b341; --bad:#e0685a; --ok:#5ac48a;
}}
@media (prefers-color-scheme: light) {{
  :root {{ --bg:#f6f7f9; --panel:#fff; --panel2:#f0f2f5; --fg:#1a1d23;
    --dim:#5a6472; --line:#e2e5ea; --accent:#1f6fb8; }}
}}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--fg);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }}
.wrap {{ max-width:1180px; margin:0 auto; padding:32px 20px 80px; }}
h1 {{ font-size:24px; margin:0 0 4px; }}
h2 {{ font-size:17px; margin:38px 0 12px; padding-bottom:6px;
  border-bottom:1px solid var(--line); }}
.meta {{ color:var(--dim); font-size:13px; }}
.pills {{ margin:14px 0 4px; display:flex; flex-wrap:wrap; gap:8px; }}
.pill {{ background:var(--panel2); border:1px solid var(--line); color:var(--dim);
  padding:3px 10px; border-radius:20px; font-size:12px; }}
.cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:12px; margin-top:18px; }}
.card {{ background:var(--panel); border:1px solid var(--line); border-radius:10px;
  padding:14px 16px; }}
.card .val {{ font-size:26px; font-weight:650; }}
.card .lbl {{ color:var(--dim); font-size:12px; margin-top:2px; }}
.card .sub {{ color:var(--dim); font-size:11px; margin-top:4px; opacity:.8; }}
.tablewrap {{ overflow-x:auto; border:1px solid var(--line); border-radius:10px;
  background:var(--panel); }}
table {{ border-collapse:collapse; width:100%; font-size:13px; }}
th,td {{ padding:7px 10px; text-align:right; white-space:nowrap;
  border-bottom:1px solid var(--line); }}
th {{ position:sticky; top:0; background:var(--panel2); color:var(--dim);
  font-weight:600; font-size:12px; text-align:right; }}
td.l, th:first-child {{ text-align:left; }}
td.dim {{ color:var(--dim); }}
td.strong {{ font-weight:700; color:var(--accent); }}
tr.warn td {{ background:color-mix(in srgb, var(--warn) 9%, transparent); }}
tbody tr:hover td {{ background:var(--panel2); }}
.note {{ color:var(--dim); font-size:13px; margin:8px 2px 12px; }}
.note.ok {{ color:var(--ok); }}
code {{ background:var(--panel2); padding:1px 5px; border-radius:4px; font-size:12px; }}
.legend {{ color:var(--dim); font-size:12px; margin:6px 2px 0; }}
th.has-tip {{ cursor:help; }}
th.has-tip::after {{ content:"\\00a0\\24d8"; color:var(--accent); font-size:10px;
  opacity:.75; font-weight:400; }}
#tip {{ position:absolute; z-index:50; display:none; max-width:340px;
  background:var(--panel); color:var(--fg); border:1px solid var(--line);
  border-radius:8px; padding:11px 13px; font-size:12px; line-height:1.5;
  text-align:left; box-shadow:0 10px 30px rgba(0,0,0,.4); pointer-events:none; }}
#tip b {{ color:var(--accent); font-size:13px; }}
#tip div {{ margin-top:5px; }}
#tip .f {{ margin-top:7px; font-family:ui-monospace,Menlo,Consolas,monospace;
  font-size:11.5px; color:var(--fg); }}
#tip .p {{ margin-top:7px; color:var(--dim); font-size:11px;
  border-top:1px solid var(--line); padding-top:6px; }}
</style></head>
<body><div class="wrap">
<h1>meta-mage &mdash; marker gene report</h1>
<div class="meta">Generated {esc(now)} &middot; {n_species:,} species across
{n_genera:,} genera &middot; {n_total:,} genomes</div>
<div class="pills">{param_items}</div>

<div class="cards">{cards}</div>

<h2>Per-species marker funnel</h2>
<p class="legend">Each stage removes gene families from the previous one:
<b>pangenome</b> &rarr; <b>core</b> (drop: below {p['core_prevalence']:.0%} in-species
prevalence) &rarr; <b>specific</b> (drop: shared with other clades / below marker
prevalence) &rarr; <b>selected</b> (drop: over the {p['max_per_clade']}/clade cap)
&rarr; <b>final</b> (drop: cross-map QC). Rows in amber ended with zero markers.</p>
{trunc_note}
<div class="tablewrap"><table><thead><tr>
<th>Species</th><th>Genus</th><th data-tip="genomes">Genomes</th>
<th data-tip="pangenome">Pangenome</th><th data-tip="core">Core</th>
<th data-tip="not_core">Not core</th><th data-tip="specific">Specific</th>
<th data-tip="not_specific">Not specific</th><th data-tip="capped">Capped</th>
<th data-tip="qc_dropped">QC dropped</th><th data-tip="final">Final</th></tr></thead>
<tbody>{''.join(body_rows)}</tbody></table></div>

<h2>Markers by rank</h2>
<div class="tablewrap"><table><thead><tr><th>Rank</th>
<th data-tip="selected">Selected</th><th data-tip="qc_dropped">QC dropped</th>
<th data-tip="final">Final</th></tr></thead>
<tbody>{''.join(rank_body)}</tbody></table></div>

<h2>Genes removed by QC (nucleotide cross-map guard)</h2>
{qc_section}

<h2>Removed species &mdash; per-genus breakdown</h2>
<p class="legend">Two ways a species leaves the DB: removed before indexing for
too few genomes, or kept but yielding zero clade-specific markers.</p>
{genus_table}

<h2>Species removed by the genome-count prefilter</h2>
{pre_section}
</div>
<script>
const TIPS = {tips_json};
(function(){{
  const tip = document.createElement('div');
  tip.id = 'tip';
  document.body.appendChild(tip);
  document.querySelectorAll('th[data-tip]').forEach(function(th){{
    const key = th.getAttribute('data-tip');
    if(!TIPS[key]) return;
    th.classList.add('has-tip');
    th.addEventListener('mouseenter', function(){{
      tip.innerHTML = TIPS[key];
      tip.style.display = 'block';
      const r = th.getBoundingClientRect();
      const w = tip.offsetWidth;
      let left = r.left + window.scrollX;
      const maxLeft = window.scrollX + document.documentElement.clientWidth - w - 12;
      if(left > maxLeft) left = maxLeft;
      if(left < window.scrollX + 8) left = window.scrollX + 8;
      tip.style.left = left + 'px';
      tip.style.top = (r.bottom + window.scrollY + 6) + 'px';
    }});
    th.addEventListener('mouseleave', function(){{ tip.style.display = 'none'; }});
  }});
}})();
</script>
</body></html>
"""


if __name__ == "__main__":
    main()
