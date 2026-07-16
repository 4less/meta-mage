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
    """Return (species_size, species_score_sum, species_genus, n_total).

    species_size[species]      = # genomes of that species (kept after prefilter)
    species_score_sum[species] = sum of those genomes' completeness (== size when
                                 no metadata was joined)
    species_genus[species]     = the genus the species belongs to
    n_total                    = total kept genomes (out-clade denominator base)
    """
    species_size = defaultdict(int)
    species_score_sum = defaultdict(float)
    species_genus = {}
    n_total = 0
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        has_comp = "completeness" in col
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sp = f[col["species"]]
            n_total += 1
            if sp == "NA":
                continue
            species_size[sp] += 1
            species_score_sum[sp] += float(f[col["completeness"]]) if has_comp else 1.0
            species_genus[sp] = f[col["genus"]]
    return species_size, species_score_sum, species_genus, n_total


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
        col = {n: i for i, n in enumerate(header.rstrip("\n").split("\t"))}
        rec_i = col.get("recovered")
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
                # A marker that cross-mapped but was kept because masking left a
                # long clean window: pass=1 AND recovered=1.
                "recovered": rec_i is not None and rec_i < len(f)
                and f[rec_i] == "1",
            }
    return verdicts


def load_merge_gain(path):
    """merge_gain.tsv rows (list of dicts) or [] when absent -- the assessment of
    whether merging a low-marker species with its neighbours lifts its markers."""
    if not path or os.path.basename(path) == "NO_FILE" or not os.path.exists(path):
        return []
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == len(header):
                rows.append(dict(zip(header, f)))
    return rows


def load_nakedness(path):
    """nakedness.tsv rows (list of dicts) or [] when absent. Classifies each
    guard-dropped marker as 'naked' (>=1 off-target clade with no competing
    marker -> truly steals reads) or 'contested' (every off-target clade markers
    the region too -> the drop was conservative)."""
    if not path or os.path.basename(path) == "NO_FILE" or not os.path.exists(path):
        return []
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == len(header):
                rows.append(dict(zip(header, f)))
    return rows


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
def build_species_funnel(counts_path, species_size, species_score_sum, n_total, p):
    """Stream counts.tsv (species rows) into a per-species funnel.

    Returns species -> counters:
        pangenome, core, non_core, specific, non_specific, below_min_in,
        clade_too_small (0/pangenome), plus the set of specific candidate keys.

    in_prevalence is completeness-weighted to match SCORE: (size - (score_sum -
    in_score)) / size, using per-species score_sum and the per-cluster in_score
    column (both default to the unweighted counts when no metadata was joined).
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
            in_score = float(f[5]) if len(f) > 5 else float(in_count)
            score_sum = species_score_sum.get(clade, float(size))
            fn = funnel[clade]
            fn["pangenome"] += 1
            in_prev = min(1.0, (size - (score_sum - in_score)) / size)
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


# Client-side boxplot renderer (kept as a plain string so its { } braces stay
# literal -- the surrounding page template is an f-string). One SVG panel per
# stage; each marker is one box (identity to every other marker in the set).
_ANI_JS = r"""
(function(){
  const STAGES = ANI_DATA.stages || [];
  const H=190, PADL=42, PADR=12, PADT=20, PADB=30;
  function yScale(v){ return PADT + (1-v)*(H-PADT-PADB); }
  function panel(boxes, title){
    const step = Math.max(6, Math.min(16, 940/Math.max(1, boxes.length)));
    const W = PADL + PADR + boxes.length*step;
    let g = '';
    [0,0.25,0.5,0.75,1].forEach(function(v){
      const y=yScale(v);
      g += '<line x1="'+PADL+'" y1="'+y+'" x2="'+W+'" y2="'+y+'" stroke="var(--line)" stroke-width="1"/>';
      g += '<text x="'+(PADL-6)+'" y="'+(y+3)+'" text-anchor="end" font-size="10" fill="var(--dim)">'+(v*100).toFixed(0)+'%</text>';
    });
    boxes.forEach(function(b,i){
      const x = PADL + i*step + step/2;
      const bw = Math.max(2, step*0.6);
      const t='<title>'+b.c+'  median '+(b.med*100).toFixed(1)+'%  max '+(b.hi*100).toFixed(1)+'%  (n='+b.n+')</title>';
      g += '<line x1="'+x+'" y1="'+yScale(b.hi)+'" x2="'+x+'" y2="'+yScale(b.lo)+'" stroke="var(--dim)" stroke-width="1"/>';
      g += '<rect x="'+(x-bw/2)+'" y="'+yScale(b.q3)+'" width="'+bw+'" height="'+Math.max(1,(yScale(b.q1)-yScale(b.q3)))+'" fill="var(--accent)" fill-opacity="0.35" stroke="var(--accent)" stroke-width="1">'+t+'</rect>';
      g += '<line x1="'+(x-bw/2)+'" y1="'+yScale(b.med)+'" x2="'+(x+bw/2)+'" y2="'+yScale(b.med)+'" stroke="var(--accent)" stroke-width="1.5"/>';
    });
    return '<div class="anipanel"><div class="anititle">'+title+' <span class="dim">(n='+boxes.length+' markers)</span></div>'
      + '<div class="aniscroll"><svg width="'+W+'" height="'+H+'">'+g+'</svg></div></div>';
  }
  function render(sp){
    const d = ANI_DATA.species[sp] || {};
    let out='';
    STAGES.forEach(function(st){ out += panel(d[st]||[], st); });
    document.getElementById('aniPanels').innerHTML = out;
  }
  const sel = document.getElementById('aniSp');
  if(sel){ sel.addEventListener('change', function(){ render(sel.value); });
    if(sel.options.length){ render(sel.value); } }
})();
"""


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
    ap.add_argument("--nakedness", help="nakedness.tsv (or NO_FILE) -- naked vs "
                    "contested classification of cross-map-dropped markers")
    ap.add_argument("--marker_ani", help="marker_ani.json (or NO_FILE) -- per-"
                    "species pairwise marker ANI at each filtering stage")
    ap.add_argument("--merge_gain", help="merge_gain.tsv (or NO_FILE) -- merge "
                    "assessment for low-marker species")
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

    species_size, species_score_sum, species_genus, n_total = load_manifest(args.manifest)
    dropped = load_dropped(args.dropped)
    verdicts = load_specificity(args.specificity)
    selected = load_selected(args.markers)
    qc_ran = len(verdicts) > 0
    funnel = build_species_funnel(args.counts, species_size, species_score_sum,
                                  n_total, p)

    # ---- per-species selected / capped / final, from the actual marker tables.
    # sp_final = markers surviving QC, INCLUDING mask-recovered ones (pass=1); of
    # those, sp_rescued were kept only because masking left a long clean window.
    sp_selected = defaultdict(int)   # species -> # selected (post-cap) markers
    sp_final = defaultdict(int)      # species -> # markers surviving QC (final)
    sp_rescued = defaultdict(int)    # species -> # of those kept by mask recovery
    for (rank, clade, cluster) in selected:
        if rank != "species":
            continue
        sp_selected[clade] += 1
        v = verdicts.get((rank, clade, cluster))
        if not qc_ran or (v and v["pass"]):
            sp_final[clade] += 1
            if v and v.get("recovered"):
                sp_rescued[clade] += 1

    # ---- markers-by-rank summary (all ranks).
    rank_selected = defaultdict(int)
    rank_final = defaultdict(int)
    rank_qc_dropped = defaultdict(int)
    rank_rescued = defaultdict(int)
    for (rank, clade, cluster) in selected:
        rank_selected[rank] += 1
        v = verdicts.get((rank, clade, cluster))
        if qc_ran and v and not v["pass"]:
            rank_qc_dropped[rank] += 1
        else:
            rank_final[rank] += 1
            if v and v.get("recovered"):
                rank_rescued[rank] += 1

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
    tot_rescued = sum(sp_rescued.values())
    n_species = len(species_size)
    n_genera = len({g for g in species_genus.values() if g != "NA"})

    merge_gain = load_merge_gain(args.merge_gain)

    # ---- nakedness of guard-dropped markers, aggregated per species.
    naked_rows = load_nakedness(args.nakedness)
    sp_naked = defaultdict(lambda: {"contested": 0, "naked": 0})
    for r in naked_rows:
        if r.get("rank") != "species":
            continue
        v = r.get("verdict")
        if v in ("contested", "naked"):
            sp_naked[r["clade"]][v] += 1
    naked_totals = {
        "contested": sum(d["contested"] for d in sp_naked.values()),
        "naked": sum(d["naked"] for d in sp_naked.values()),
    }

    # ---- per-species pairwise marker ANI (report boxplots), embedded verbatim.
    marker_ani = {}
    if (args.marker_ani and os.path.basename(args.marker_ani) != "NO_FILE"
            and os.path.exists(args.marker_ani)):
        with open(args.marker_ani) as fh:
            marker_ani = json.load(fh)

    html_out = render_html(
        args, p, qc_ran, n_total, n_species, n_genera,
        species_size, species_genus, funnel, sp_selected, sp_final, sp_rescued,
        rank_selected, rank_final, rank_qc_dropped, rank_rescued,
        qc_removed, dropped, genus_species, genus_prefilter, genus_zero,
        zero_marker_species, merge_gain, sp_naked, naked_totals, marker_ani,
        totals=dict(pangenome=tot_pangenome, core=tot_core, specific=tot_specific,
                    selected=tot_selected, final=tot_final, rescued=tot_rescued),
    )
    with open(args.out, "w") as fh:
        fh.write(html_out)
    print(f"build_report: wrote {args.out} "
          f"({n_species} species, {tot_final} final markers, "
          f"{len(qc_removed)} removed by QC)")


def render_html(args, p, qc_ran, n_total, n_species, n_genera,
                species_size, species_genus, funnel, sp_selected, sp_final, sp_rescued,
                rank_selected, rank_final, rank_qc_dropped, rank_rescued,
                qc_removed, dropped, genus_species, genus_prefilter, genus_zero,
                zero_marker_species, merge_gain, sp_naked, naked_totals,
                marker_ani, totals):
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # Summary cards.
    rescued = totals.get("rescued", 0)
    final_sub = ("QC not run" if not qc_ran else
                 f"after cross-map QC, incl. {rescued:,} mask-rescued" if rescued
                 else "after cross-map QC + mask recovery")
    cards = "".join([
        card("Genomes", f"{n_total:,}"),
        card("Species", f"{n_species:,}"),
        card("Genera", f"{n_genera:,}"),
        card("Gene families", f"{totals['pangenome']:,}", "species-rank, summed"),
        card("Core genes", f"{totals['core']:,}", f">= {p['core_prevalence']:.0%} prevalence"),
        card("Specific candidates", f"{totals['specific']:,}", "pre-cap"),
        card("Selected markers", f"{totals['selected']:,}", f"cap {p['max_per_clade']}/clade"),
        card("Final markers", f"{totals['final']:,}", final_sub),
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
        resc = sp_rescued.get(sp, 0)
        capped = max(0, spec - sel)
        # QC dropped = selected minus final; final already includes mask-rescued,
        # so a marker that cross-mapped but was recovered is NOT counted here.
        qc_drop = max(0, sel - fin) if qc_ran else 0
        cls = " class='warn'" if fin == 0 else ""
        # Survivor columns are plain; loss columns carry .loss (dim).
        qc_cell = f"<td class='loss'>{qc_drop}</td>" if qc_ran else "<td class='loss'>-</td>"
        resc_cell = f"<td>{resc}</td>" if qc_ran else "<td>-</td>"
        body_rows.append(
            f"<tr{cls}><td class='l'>{esc(sp)}</td>"
            f"<td class='l dim'>{esc(species_genus.get(sp, 'NA'))}</td>"
            f"<td>{species_size.get(sp, 0)}</td>"
            f"<td>{pan}</td>"
            f"<td class='loss'>{fn.get('non_core', 0)}</td><td>{core}</td>"
            f"<td class='loss'>{fn.get('non_specific', 0) + fn.get('below_min_in', 0)}</td>"
            f"<td>{spec}</td>"
            f"<td class='loss'>{capped}</td><td>{sel}</td>"
            f"{qc_cell}{resc_cell}"
            f"<td class='strong'>{fin}</td></tr>"
        )
    trunc_note = (f"<p class='note'>Showing the first {args.max_table_rows:,} of "
                  f"{len(species_rows):,} species (sorted by genus). Full per-marker "
                  f"data is in the pipeline's counts/markers TSVs.</p>"
                  if truncated else "")

    # Merge assessment: for low-marker species, whether merging with the nearest
    # neighbours (markers recomputed after all filters, masking RE-APPLIED for the
    # merged clade) lifts the marker count, and which species it would pool.
    merge_section = ""
    if merge_gain:
        mg = sorted(merge_gain,
                    key=lambda r: (r.get("reached_threshold") != "yes",
                                   -int(r.get("delta", 0) or 0)))
        n_reach = sum(1 for r in mg if r.get("reached_threshold") == "yes")
        n_help = sum(1 for r in mg if int(r.get("delta", 0) or 0) > 0)
        mrows = []
        for r in mg[: args.max_table_rows]:
            reached = r.get("reached_threshold") == "yes"
            delta = int(r.get("delta", 0) or 0)
            added = int(r.get("n_species_added", 0) or 0)
            merged_with = " + ".join(
                (r.get("merged_clade", "") or "").split(" + ")[1:]) or "&mdash;"
            badge = ("<span class='ok'>reaches target</span>" if reached
                     else f"<span class='warn'>+{delta}</span>" if delta > 0
                     else "<span class='dim'>no gain</span>")
            cls = "" if delta > 0 else " class='dim'"
            mrows.append(
                f"<tr{cls}><td class='l'>{esc(r.get('species',''))}</td>"
                f"<td class='l dim'>{esc(r.get('genus',''))}</td>"
                f"<td>{esc(r.get('baseline_markers','0'))}</td>"
                f"<td class='strong'>{esc(r.get('merged_markers','0'))}</td>"
                f"<td>{delta:+d}</td><td>{added}</td>"
                f"<td class='l'>{merged_with}</td>"
                f"<td>{esc(r.get('trajectory',''))}</td><td>{badge}</td></tr>")
        merge_section = (
            "<h2>Merge assessment (low-marker species)</h2>"
            "<p class='legend'>For each species that ended below the marker "
            "threshold, the greedy merge probe recomputes markers <b>after all "
            "filters, re-applying masking for the merged clade</b>, absorbing the "
            "nearest same-genus neighbours until the count clears the threshold. "
            f"<b>{n_help}</b> species gain markers by merging, <b>{n_reach}</b> reach "
            "the threshold. Baseline = current final markers; Merged = after the "
            "listed species are pooled.</p>"
            "<div class='tablewrap'><table><thead><tr>"
            "<th>Species</th><th>Genus</th><th>Final markers</th>"
            "<th>Merged markers</th><th>&Delta;</th><th># merged</th>"
            "<th>Merge with</th><th>Trajectory</th><th>Verdict</th></tr></thead>"
            f"<tbody>{''.join(mrows)}</tbody></table></div>")

    # Nakedness of cross-map drops: naked (no competing off-target marker ->
    # truly steals reads) vs contested (off-target clade markers the region too
    # -> the guard dropped it conservatively; a competitive rule could keep it).
    naked_section = ""
    if sp_naked:
        tot_c = naked_totals["contested"]
        tot_n = naked_totals["naked"]
        tot = tot_c + tot_n
        nrows = sorted(
            sp_naked.items(),
            key=lambda kv: (-(kv[1]["contested"] + kv[1]["naked"]),
                            -kv[1]["contested"]),
        )
        nbody = []
        for sp, d in nrows[: args.max_table_rows]:
            drp = d["contested"] + d["naked"]
            cls = " class='warn'" if d["contested"] and d["naked"] == 0 else ""
            nbody.append(
                f"<tr{cls}><td class='l'>{esc(sp)}</td>"
                f"<td class='l dim'>{esc(species_genus.get(sp, 'NA'))}</td>"
                f"<td>{drp}</td><td class='strong'>{d['contested']}</td>"
                f"<td>{d['naked']}</td><td>{pct(d['naked'], drp)}</td></tr>"
            )
        body = "".join(nbody)
        naked_section = (
            "<h2>Cross-map drops: naked vs contested</h2>"
            "<p class='legend'>The guard drops a marker whenever its sequence "
            "occurs outside its clade. But a drop only truly costs specificity "
            "when the off-target region has <b>no competing marker</b> on the "
            "other side (<b>naked</b> &mdash; reads land only here, inflating the "
            "wrong taxon). When the off-target clade markers the same region "
            "(<b>contested</b>), reads are shared, so under a competitive "
            "assignment rule the marker could have been kept. Contested drops are "
            "where the conservative guard is over-dropping. Of "
            f"<b>{tot:,}</b> classified cross-map drops, <b>{tot_c:,}</b> are "
            f"contested (recoverable) and <b>{tot_n:,}</b> naked (correctly "
            "dropped). Amber rows lost only contested markers.</p>"
            "<div class='tablewrap'><table><thead><tr><th>Species</th><th>Genus</th>"
            "<th>Cross-map dropped</th><th>Contested (recoverable)</th>"
            "<th>Naked</th><th>% naked</th></tr></thead>"
            f"<tbody>{body}</tbody></table></div>"
        )

    # Per-species pairwise marker ANI boxplots (one panel per filtering stage;
    # each box = a marker's identity to every other marker in that stage's set).
    ani_section = ""
    if marker_ani and marker_ani.get("species"):
        sp_data = marker_ani["species"]
        # default + option order: most stage-1 markers first (most to look at).
        sp_keys = sorted(sp_data,
                         key=lambda s: -len(sp_data[s].get("specific-200", [])))
        options = "".join(
            "<option value='%s'>%s</option>" % (esc(s), esc(s)) for s in sp_keys)
        data_json = json.dumps(marker_ani).replace("</", "<\\/")
        ani_section = (
            "<h2>Pairwise marker ANI per species</h2>"
            "<p class='legend'>For a species' markers, each <b>box is one marker</b>"
            " and summarises its pairwise nucleotide identity (mash-style, "
            "k=" + str(marker_ani.get("kmer", "?")) + ") to <b>every other marker"
            "</b> in that stage's set. Low boxes = the marker is distinct; a high "
            "box or tall whisker = it is redundant with (or cross-similar to) other "
            "markers. Panels are the three filtering stages, each capped at "
            + str(marker_ani.get("cap", 200)) + " markers ordered by score: "
            "<b>specific-200</b> (top candidates, pre-guard) &rarr; "
            "<b>post-crossmap</b> (survived the guard cleanly) &rarr; "
            "<b>post-recovery</b> (final, incl. mask-rescued).</p>"
            "<div class='anictl'>Species: <select id='aniSp'>" + options
            + "</select></div><div id='aniPanels'></div>"
            "<script>\nconst ANI_DATA=" + data_json + ";\n" + _ANI_JS
            + "\n</script>"
        )

    # Markers-by-rank table.
    rank_body = []
    for rank in RANKS:
        if rank_selected.get(rank, 0) == 0:
            continue
        rank_body.append(
            f"<tr><td class='l'>{esc(rank)}</td>"
            f"<td>{rank_selected.get(rank, 0)}</td>"
            f"<td>{rank_qc_dropped.get(rank, 0) if qc_ran else '-'}</td>"
            f"<td>{rank_rescued.get(rank, 0) if qc_ran else '-'}</td>"
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
td.loss, th.loss {{ color:var(--dim); font-weight:400;
  background:color-mix(in srgb, var(--line) 30%, transparent); }}
.anictl {{ margin:10px 2px; color:var(--dim); font-size:13px; }}
.anictl select {{ background:var(--panel2); color:var(--fg);
  border:1px solid var(--line); border-radius:6px; padding:4px 8px; font-size:13px; }}
.anipanel {{ margin:14px 0; }}
.anititle {{ font-size:13px; color:var(--fg); margin-bottom:2px; }}
.aniscroll {{ overflow-x:auto; border:1px solid var(--line);
  border-radius:8px; background:var(--panel); padding:4px 0; }}
td.strong {{ font-weight:700; color:var(--accent); }}
span.ok {{ color:var(--ok); font-weight:600; }}
span.warn {{ color:var(--warn); font-weight:600; }}
span.dim {{ color:var(--dim); }}
tr.dim td {{ color:var(--dim); }}
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
<p class="legend">Bold columns are the <b>surviving marker count at each stage</b>;
dim columns are what was <b>lost</b> to reach the next survivor count. Reading left
to right: <b>Pan</b> &rarr; (&minus;<b>&not;Core</b>) &rarr; <b>Core</b> &rarr;
(&minus;<b>&not;Spec</b>) &rarr; <b>Spec</b> &rarr; (&minus;<b>Cap</b>, the
{p['max_per_clade']}/clade cap) &rarr; <b>Sel</b> &rarr; (&minus;<b>xMap&#10007;</b>
cross-map drops, &plus;<b>Resc</b> mask-rescued) &rarr; <b>Fin</b>. So
Fin&nbsp;=&nbsp;Sel&nbsp;&minus;&nbsp;xMap&#10007;, and Resc is the rescued subset
already inside Fin. Hover any header for its formula. Rows in amber ended with zero
markers.</p>
{trunc_note}
<div class="tablewrap"><table><thead><tr>
<th>Species</th><th>Genus</th><th data-tip="genomes">N</th>
<th data-tip="pangenome">Pan</th>
<th data-tip="not_core" class="loss">&not;Core</th><th data-tip="core">Core</th>
<th data-tip="not_specific" class="loss">&not;Spec</th><th data-tip="specific">Spec</th>
<th data-tip="capped" class="loss">Cap</th><th data-tip="selected">Sel</th>
<th data-tip="qc_dropped" class="loss">xMap&#10007;</th>
<th data-tip="rescued">Resc</th>
<th data-tip="final">Fin</th></tr></thead>
<tbody>{''.join(body_rows)}</tbody></table></div>

{ani_section}

{merge_section}

<h2>Markers by rank</h2>
<div class="tablewrap"><table><thead><tr><th>Rank</th>
<th data-tip="selected">Selected</th><th data-tip="qc_dropped">QC dropped</th>
<th data-tip="rescued">Mask-rescued</th>
<th data-tip="final">Final</th></tr></thead>
<tbody>{''.join(rank_body)}</tbody></table></div>

<h2>Genes removed by QC (nucleotide cross-map guard)</h2>
{qc_section}

{naked_section}

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
