from __future__ import annotations

from pathlib import Path

import pandas as pd

from aadr_parser import load_aadr_table


def normalize_y_label(label: str) -> str:
    if label is None:
        return ""
    return (
        str(label)
        .strip()
        .replace("?", "")
        .replace("*", "")
    )


def export_y_aliases(aadr_path: str, output_path: str) -> None:
    df = load_aadr_table(aadr_path)

    required_cols = ["y_haplogroup", "y_haplogroup_isogg"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column in AADR data: {col}")

    ymap = df[["y_haplogroup", "y_haplogroup_isogg"]].copy()

    ymap["y_haplogroup"] = ymap["y_haplogroup"].apply(normalize_y_label)
    ymap["y_haplogroup_isogg"] = ymap["y_haplogroup_isogg"].apply(normalize_y_label)

    ymap = ymap[
        (ymap["y_haplogroup"] != "") &
        (ymap["y_haplogroup"].notna()) &
        (ymap["y_haplogroup_isogg"] != "") &
        (ymap["y_haplogroup_isogg"].notna())
    ].drop_duplicates()

    # Export as an alias table format
    alias_df = pd.DataFrame({
        "raw_label": ymap["y_haplogroup"],
        "aadr_isogg": ymap["y_haplogroup_isogg"],
        "tree_label": ymap["y_haplogroup_isogg"],
    })

    alias_df = alias_df.drop_duplicates().sort_values(
        ["raw_label", "aadr_isogg"]
    ).reset_index(drop=True)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    alias_df.to_csv(output_path, sep="\t", index=False)

    print(f"Exported {len(alias_df)} Y alias rows to: {output_path}")


if __name__ == "__main__":
    aadr_file = "../data/aadr/AADR Annotations 2025.xlsx"
    out_file = "../data/y_aliases.tsv"
    export_y_aliases(aadr_file, out_file)
