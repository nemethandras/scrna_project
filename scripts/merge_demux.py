"""
Merge Vireo assignments with DB donor matching and binomial scorer output
into a single final assignments table.

For each cell barcode:
  - If Vireo's donor maps to a DB reference (concordance >= threshold)
    → cell_line name from DB, source = db_match
  - If Vireo's donor has no DB match
    → labelled "vireo:<donor_id>", source = vireo_only
  - If Vireo flagged as doublet or unassigned
    → kept as-is

Optionally joins the binomial scorer (demux.py) output for comparison.

Usage:
    python scripts/merge_demux.py \
        --vireo-dir   results/demux/pool_ctr_001/vireo \
        --matches     results/demux/pool_ctr_001/donor_matches.tsv \
        --output      results/demux/pool_ctr_001/final_assignments.tsv \
        [--scorer     results/demux/pool_ctr_001/assignments.tsv] \
        [--min-concordance 0.9]
"""

import argparse
import csv
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Merge Vireo + DB-match + binomial scorer into final assignments"
    )
    p.add_argument("--vireo-dir",       required=True,
                   help="Vireo output directory (must contain donor_ids.tsv)")
    p.add_argument("--matches",         required=True,
                   help="donor_matches.tsv from match_vireo.py")
    p.add_argument("--output",          required=True,
                   help="output TSV path for final merged assignments")
    p.add_argument("--scorer",          default=None,
                   help="assignments.tsv from demux.py (optional, adds comparison column)")
    p.add_argument("--min-concordance", type=float, default=0.65, metavar="F",
                   help="ignored — kept for CLI compatibility; acceptance is driven by the "
                        "confidence column written by match_vireo.py (high/low = accept, "
                        "no_match/no_data = reject)")
    return p.parse_args()


def load_donor_matches(matches_path, min_concordance):
    """
    Returns dict: vireo_donor → (cell_line, concordance, n_positions, confidence) or None.
    Uses the cell_line column (friendly name) when available, else assigned_line (run_id).
    """
    mapping = {}
    with open(matches_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            donor      = row["vireo_donor"]
            run_id     = row["assigned_line"]
            cell_line  = row.get("cell_line", "").strip() or run_id
            conc       = float(row["concordance"]) if row["concordance"] else 0.0
            n          = int(row["n_positions"])   if row["n_positions"]  else 0
            confidence = row.get("confidence", "").strip()
            # Accept high and low confidence — match_vireo.py already applied
            # the threshold; no_match / no_data means genuinely unassigned.
            if run_id != "no_match" and n > 0 and confidence in ("high", "low"):
                mapping[donor] = (cell_line, conc, n, confidence)
            else:
                mapping[donor] = None
    return mapping


def load_vireo_assignments(vireo_dir):
    """
    Returns dict: barcode → {donor_id, prob_max, prob_doublet, n_vars}
    """
    path = Path(vireo_dir) / "donor_ids.tsv"
    if not path.exists():
        raise SystemExit(f"\nERROR: {path} not found — did Vireo complete successfully?\n")
    cells = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cells[row["cell"]] = row
    return cells


def load_scorer_assignments(scorer_path):
    """Returns dict: barcode → row from demux.py assignments.tsv"""
    if not scorer_path or not Path(scorer_path).exists():
        return {}
    result = {}
    with open(scorer_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            result[row["barcode"]] = row
    return result


def main():
    args = parse_args()

    print(f"Loading Vireo assignments from: {args.vireo_dir}")
    vireo = load_vireo_assignments(args.vireo_dir)
    print(f"  {len(vireo):,} cells")

    print(f"Loading donor→cell_line matches from: {args.matches}")
    donor_map = load_donor_matches(args.matches, args.min_concordance)
    print(f"  {sum(v is not None for v in donor_map.values())}/{len(donor_map)} donors matched to DB "
          f"(concordance >= {args.min_concordance})")
    for donor, match in sorted(donor_map.items()):
        if match:
            line, conc, n, conf = match
            print(f"    {donor:<12} → {line}  (concordance={conc:.3f}, n={n:,}, confidence={conf})")
        else:
            print(f"    {donor:<12} → no_match")

    scorer = load_scorer_assignments(args.scorer)
    has_scorer = bool(scorer)

    # Build final assignments
    counts = {"db_match": 0, "vireo_only": 0, "doublet": 0, "unassigned": 0}

    with open(args.output, "w") as f:
        header = [
            "barcode", "cell_line", "source",
            "vireo_donor", "vireo_prob_max", "vireo_prob_doublet",
            "concordance", "n_positions", "match_confidence",
        ]
        if has_scorer:
            header += ["scorer_assignment", "scorer_status", "agreement"]
        f.write("\t".join(header) + "\n")

        for barcode, vrow in vireo.items():
            donor    = vrow["donor_id"]
            prob_max = vrow.get("prob_max", "")
            prob_dbl = vrow.get("prob_doublet", "")

            if donor == "doublet":
                # Resolve constituent donors to cell line names where possible
                best_doublet = vrow.get("best_doublet", "")
                if best_doublet:
                    parts = [p.strip() for p in best_doublet.split(",")]
                    resolved = []
                    for part in parts:
                        m = donor_map.get(part)
                        resolved.append(m[0] if m else f"vireo:{part}")
                    cell_line = "doublet:" + "+".join(resolved)
                else:
                    cell_line = "doublet"
                source      = "doublet"
                concordance = ""
                n_pos       = ""
                confidence  = ""
                counts["doublet"] += 1
            elif donor == "unassigned":
                cell_line   = "unassigned"
                source      = "unassigned"
                concordance = ""
                n_pos       = ""
                confidence  = ""
                counts["unassigned"] += 1
            else:
                match = donor_map.get(donor)
                if match:
                    cell_line, conc, n, conf = match
                    concordance = f"{conc:.4f}"
                    n_pos       = str(n)
                    confidence  = conf
                    source      = "db_match"
                    counts["db_match"] += 1
                else:
                    cell_line   = f"vireo:{donor}"
                    concordance = ""
                    n_pos       = ""
                    confidence  = ""
                    source      = "vireo_only"
                    counts["vireo_only"] += 1

            row = [barcode, cell_line, source, donor, prob_max, prob_dbl, concordance, n_pos, confidence]

            if has_scorer:
                srow       = scorer.get(barcode, {})
                scorer_asgn  = srow.get("assigned_line", "")
                scorer_status = srow.get("status", "")
                # Agreement: both assign the same cell line (ignoring source prefix)
                agree = str(scorer_asgn == cell_line).upper() if scorer_asgn else ""
                row += [scorer_asgn, scorer_status, agree]

            f.write("\t".join(row) + "\n")

    total = len(vireo)
    print(f"\nFinal assignments ({total:,} cells):")
    for status, count in counts.items():
        print(f"  {status:<12}  {count:>6,}  ({100*count/total:.1f}%)")
    print(f"\nOutput: {args.output}")


if __name__ == "__main__":
    main()
