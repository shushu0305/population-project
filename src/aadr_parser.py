from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import pandas as pd


def load_aadr_table(
    file_path: str | Path,
    sep: str = "\t",
    sheet_name: Optional[Union[str, int]] = 0,
) -> pd.DataFrame:
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"AADR file not found: {file_path}")

    if file_path.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(file_path, sheet_name=sheet_name)
    else:
        df = pd.read_csv(file_path, sep=sep, low_memory=False)

    df.columns = [str(col).strip() for col in df.columns]

    rename_map = {
        "Master ID": "sample_id",
        "Genetic ID": "genetic_id",

        "mtDNA haplogroup": "mt_haplogroup",
        "mtDNA haplogroup if >2x or published": "mt_haplogroup",

        "Y haplogroup (manual curation in terminal mutation format)": "y_haplogroup",
        "Y haplogroup (manual curation in ISOGG format)": "y_haplogroup_isogg",

        "Date mean in BP": "date_mean_bp",
        "Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]": "date_mean_bp",

        "Political Entity": "political_entity",
        "Locality": "locality",
        "Publication": "publication",
        "Group ID": "group_id",
        "Libraries": "libraries",
    }

    df = df.rename(columns={c: rename_map.get(c, c) for c in df.columns})

    needed_cols = [
        "sample_id",
        "mt_haplogroup",
        "y_haplogroup",
        "date_mean_bp",
        "political_entity",
        "locality",
    ]

    optional_cols = [
        "genetic_id",
        "y_haplogroup_isogg",
        "publication",
        "group_id",
        "libraries",
    ]

    for col in needed_cols + optional_cols:
        if col not in df.columns:
            df[col] = None

    # Clean the string field
    string_cols = [
        "sample_id",
        "mt_haplogroup",
        "y_haplogroup",
        "y_haplogroup_isogg",
        "political_entity",
        "locality",
        "genetic_id",
        "publication",
        "group_id",
        "libraries",
    ]

    for col in string_cols:
        df[col] = (
            df[col]
            .astype(str)
            .str.strip()
            .replace({"..": None, "": None, "nan": None, "n/a": None})
        )

    # Remove the common uncertain symbols in mt/Y
    for hg_col in ["mt_haplogroup", "y_haplogroup", "y_haplogroup_isogg"]:
        df[hg_col] = (
            df[hg_col]
            .astype(str)
            .str.strip()
            .str.replace("?", "", regex=False)
            .str.replace("*", "", regex=False)
            .replace({"None": None, "nan": None})
        )

    # convert Date to value
    df["date_mean_bp"] = pd.to_numeric(df["date_mean_bp"], errors="coerce")

    return df[needed_cols + optional_cols].copy()