"""
main.py
-------
Command-line interface (CLI) entry point for Haplogroup Discover.

Usage examples
~~~~~~~~~~~~~~
# mtDNA query, print results to terminal
python main.py --system mt --target U5b2c \\
    --aadr "../data/aadr/AADR Annotations 2025.xlsx"

# Y-DNA query with VIP matching, save output tables and report
python main.py --system y --target R-M269 \\
    --aadr "../data/aadr/AADR Annotations 2025.xlsx" \\
    --vip ../data/vip/VIPHaplogroups.xlsx \\
    --early-n 10 --save

Run with -h / --help for a full argument reference.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

from aadr_parser import load_aadr_table
from analysis_origin import HaplogroupAnalyzer as OriginHaplogroupAnalyzer
from analysis_vip import HaplogroupAnalyzer as VipHaplogroupAnalyzer
from tree_parser import PhyloTree
from vip_parser import load_vip_table
from config import DATA_CONFIG


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze ancient haplogroup occurrence, origin candidate, and VIP matches."
    )

    parser.add_argument(
        "--system",
        required=True,
        choices=["mt", "y"],
        help="Haplogroup system: mt or y",
    )
    parser.add_argument(
        "--target",
        required=True,
        help="Target haplogroup, e.g. U5b2c or R-M269",
    )
    parser.add_argument(
        "--aadr",
        required=True,
        help="Path to AADR annotation file",
    )
    parser.add_argument(
        "--vip",
        required=False,
        default=None,
        help="Path to VIP Excel file",
    )
    parser.add_argument(
        "--vip-sheet",
        default=None,
        help="VIP Excel sheet name/index. Default: auto-select based on --system, or load all sheets when unclear.",
    )
    parser.add_argument(
        "--mt-tree",
        default=DATA_CONFIG["MT_TREE_PATH"],
        help="Path to mtDNA tree file",
    )
    parser.add_argument(
        "--y-tree",
        default=DATA_CONFIG["Y_TREE_PATH"],
        help="Path to Y tree file",
    )
    parser.add_argument(
        "--early-n",
        type=int,
        default=5,
        help="Number of earliest samples used for early-country summary",
    )
    parser.add_argument(
        "--output-dir",
        default="output",
        help="Directory for saving results",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save result tables to files",
    )

    return parser.parse_args()


def ensure_output_dirs(base_dir: Path) -> dict[str, Path]:
    tables_dir = base_dir / "tables"
    reports_dir = base_dir / "reports"

    tables_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    return {
        "base": base_dir,
        "tables": tables_dir,
        "reports": reports_dir,
    }


def load_trees(mt_tree_path: str, y_tree_path: str) -> tuple[PhyloTree, PhyloTree]:
    mt_tree = PhyloTree(mt_tree_path)
    mt_tree.load()

    y_tree = PhyloTree(y_tree_path)
    y_tree.load()

    return mt_tree, y_tree


def combine_results(origin_result, vip_result):
    """Keep origin/early-sample outputs from origin_result and VIP from vip_result."""
    return SimpleNamespace(
        target=origin_result.target,
        system=origin_result.system,
        resolved_target_for_tree=getattr(origin_result, "resolved_target_for_tree", None),
        upstream_label=origin_result.upstream_label,
        downstream_labels=getattr(origin_result, "downstream_labels", []),
        included_labels=origin_result.included_labels,
        matched_samples=origin_result.matched_samples,
        oldest_sample=origin_result.oldest_sample,
        early_samples=origin_result.early_samples,
        country_summary=origin_result.country_summary,
        candidate_origin_country=origin_result.candidate_origin_country,
        vip_matches=vip_result.vip_matches,
    )


def save_result_tables(result, output_dirs: dict[str, Path]) -> None:
    safe_target = result.target.replace("/", "_").replace(" ", "_")
    prefix = f"{result.system}_{safe_target}"

    result.matched_samples.to_csv(
        output_dirs["tables"] / f"{prefix}_matched_samples.csv",
        index=False,
    )
    result.oldest_sample.to_csv(
        output_dirs["tables"] / f"{prefix}_oldest_sample.csv",
        index=False,
    )
    result.early_samples.to_csv(
        output_dirs["tables"] / f"{prefix}_early_samples.csv",
        index=False,
    )
    result.country_summary.to_csv(
        output_dirs["tables"] / f"{prefix}_country_summary.csv",
        index=False,
    )
    result.vip_matches.to_csv(
        output_dirs["tables"] / f"{prefix}_vip_matches.csv",
        index=False,
    )

    report_path = output_dirs["reports"] / f"{prefix}_summary.txt"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(build_summary_text(result))


def build_summary_text(result) -> str:
    lines = []
    lines.append(f"Target: {result.target}")
    lines.append(f"System: {result.system}")
    lines.append(f"Direct upstream label: {result.upstream_label}")
    lines.append("")

    lines.append("Included labels:")
    for label in result.included_labels[:50]:
        lines.append(f"  - {label}")
    if len(result.included_labels) > 50:
        lines.append(f"  ... ({len(result.included_labels)} total)")
    lines.append("")

    lines.append("Oldest sample:")
    if result.oldest_sample.empty:
        lines.append("  No matched samples found.")
    else:
        row = result.oldest_sample.iloc[0]
        hg_value = row.get("mt_haplogroup") if result.system == "mt" else row.get("y_haplogroup")
        lines.append(f"  sample_id: {row.get('sample_id')}")
        lines.append(f"  haplogroup: {hg_value}")
        lines.append(f"  date_mean_bp: {row.get('date_mean_bp')}")
        lines.append(f"  political_entity: {row.get('political_entity')}")
        lines.append(f"  locality: {row.get('locality')}")
    lines.append("")

    lines.append("Candidate origin country:")
    lines.append(f"  {result.candidate_origin_country}")
    lines.append("")

    lines.append("VIP matches:")
    if result.vip_matches.empty:
        lines.append("  No VIP matches found.")
    else:
        for _, row in result.vip_matches.iterrows():
            lines.append(
                f"  - {row['vip_name']} | {row['vip_haplogroup']} | {row['relation']}"
            )

    return "\n".join(lines)


def print_result(result) -> None:
    print("=== Direct upstream label ===")
    print(result.upstream_label)

    print("\n=== Included labels ===")
    print(result.included_labels[:30])
    if len(result.included_labels) > 30:
        print(f"... ({len(result.included_labels)} total)")

    print("\n=== Matched samples ===")
    print(result.matched_samples)

    print("\n=== Oldest sample ===")
    print(result.oldest_sample)

    print("\n=== Early samples ===")
    print(result.early_samples)

    print("\n=== Country summary ===")
    print(result.country_summary)

    print("\n=== Candidate origin country ===")
    print(result.candidate_origin_country)

    print("\n=== VIP matches ===")
    print(result.vip_matches)


def main() -> None:
    args = parse_args()

    mt_tree, y_tree = load_trees(args.mt_tree, args.y_tree)

    aadr_df = load_aadr_table(args.aadr, sep="\t")

    vip_df = None
    if args.vip:
        vip_sheet = args.vip_sheet
        if isinstance(vip_sheet, str) and vip_sheet.isdigit():
            vip_sheet = int(vip_sheet)

        if vip_sheet is None:
            vip_sheet = "Y" if args.system == "y" else "mtDNA"

        vip_df = load_vip_table(args.vip, sheet_name=vip_sheet)

    origin_analyzer = OriginHaplogroupAnalyzer(
        aadr_df=aadr_df,
        mt_tree=mt_tree,
        y_tree=y_tree,
        vip_df=vip_df,
    )
    vip_analyzer = VipHaplogroupAnalyzer(
        aadr_df=aadr_df,
        mt_tree=mt_tree,
        y_tree=y_tree,
        vip_df=vip_df,
    )

    origin_result = origin_analyzer.analyze(
        target=args.target,
        system=args.system,
        early_n=args.early_n,
    )
    vip_result = vip_analyzer.analyze(
        target=args.target,
        system=args.system,
        early_n=args.early_n,
    )

    result = combine_results(origin_result, vip_result)

    print_result(result)

    if args.save:
        output_dirs = ensure_output_dirs(Path(args.output_dir))
        save_result_tables(result, output_dirs)
        print(f"\nSaved results to: {output_dirs['base']}")


if __name__ == "__main__":
    main()
