#!/usr/bin/env python3
"""
Pick the relaxed-tier candidates worth guarding for the adaptive rescue.

From the relaxed scored table (markers scored at the low in-prevalence floor),
keep only the rows that are:
  * for a FLAGGED species (still below the marker threshold after mask recovery,
    and whose merge probe does not reach the threshold), AND
  * below the normal in-prevalence start (in_prev < min_in_start) -- i.e. the NEW
    lower-prevalence candidates that the 0.8 run never selected.

These are the extra candidates fed to EMIT_REPS -> CROSSMAP -> clean guard.
"""
import argparse


def load_flagged(low_marker, merge_gain):
    """Species that are flagged (below threshold) AND whose merge does not reach
    the threshold -- the ones the adaptive relaxation should try to rescue."""
    flagged = set()
    with open(low_marker) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ci["flagged"]] == "yes":
                flagged.add(f[ci["species"]])
    # Drop species a merge would rescue (merge is preferred; don't also relax them).
    if merge_gain:
        try:
            fh = open(merge_gain)
        except OSError:
            fh = None
        if fh:
            with fh:
                header = fh.readline().rstrip("\n").split("\t")
                ci = {n: i for i, n in enumerate(header)}
                si, ri = ci.get("species"), ci.get("reached_threshold")
                if si is not None and ri is not None:
                    for line in fh:
                        f = line.rstrip("\n").split("\t")
                        if len(f) > max(si, ri) and f[ri] == "yes":
                            flagged.discard(f[si])
    return flagged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--relaxed", required=True, help="markers_relaxed.tsv (scored at floor)")
    ap.add_argument("--low_marker", required=True, help="low_marker_species.tsv")
    ap.add_argument("--merge_gain", help="merge_gain.tsv (or NO_FILE)")
    ap.add_argument("--rank", default="species")
    ap.add_argument("--min_in_start", type=float, required=True,
                    help="normal in-prevalence (rows >= this were already selected)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    mg = args.merge_gain
    if mg and mg.endswith("NO_FILE"):
        mg = None
    flagged = load_flagged(args.low_marker, mg)

    kept = 0
    with open(args.relaxed) as fh, open(args.out, "w") as out:
        header = fh.readline()
        out.write(header)
        cols = header.rstrip("\n").split("\t")
        ip_i = cols.index("in_prevalence")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != args.rank or f[1] not in flagged:
                continue
            try:
                if float(f[ip_i]) >= args.min_in_start:
                    continue  # already covered by the normal 0.8 run
            except (ValueError, IndexError):
                continue
            out.write(line if line.endswith("\n") else line + "\n")
            kept += 1
    print(f"filter_relax: {len(flagged)} flagged species, {kept} extra candidates")


if __name__ == "__main__":
    main()
