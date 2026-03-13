from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import pandas as pd


VIP_Y_SHEET_NAMES = {"y", "ychrom", "y-chrom", "y_chrom", "y chromosome", "y-chromosome"}
VIP_MT_SHEET_NAMES = {"mtdna", "mt", "mtdna "}


def _normalize_vip_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(col).strip() for col in df.columns]

    rename_map = {
        "Individual": "vip_name",
        "Name": "vip_name",
        "VIP": "vip_name",
        "Category": "category",
        "mtDNA": "haplogroup",
        "Haplogroup": "haplogroup",
        "HG": "haplogroup",
        "System": "system",
        "system": "system",
        "Source": "source",
        "source": "source",
        "Note": "note",
        "note": "note",
    }

    df = df.rename(columns={c: rename_map.get(c, c) for c in df.columns})

    if "vip_name" not in df.columns:
        raise ValueError("VIP file missing required column: vip_name")
    if "haplogroup" not in df.columns:
        raise ValueError("VIP file missing required column: haplogroup")

    for col in ["system", "source", "note", "category"]:
        if col not in df.columns:
            df[col] = None

    df["vip_name"] = df["vip_name"].astype(str).str.strip()
    df["haplogroup"] = df["haplogroup"].astype(str).str.strip()

    return df


def _infer_system_from_sheet_name(sheet_name: Optional[Union[str, int]]) -> Optional[str]:
    if not isinstance(sheet_name, str):
        return None
    s = sheet_name.strip().lower()
    if s in VIP_MT_SHEET_NAMES:
        return "mt"
    if s in VIP_Y_SHEET_NAMES:
        return "y"
    return None


def _infer_system_from_hg(hg: str) -> str:
    hg = str(hg or "").strip()
    if not hg:
        return ""

    if "-" in hg:
        return "y"

    y_only_prefixes = ("A0", "A1", "BT", "CT", "DE", "CF")
    if any(hg.startswith(p) for p in y_only_prefixes):
        return "y"

    # Y ISOGG-style labels such as R1b1a1b1a1a2b1, E1b1a, J2a1, etc.
    # Heuristic: one or two leading uppercase letters followed by digits,
    # usually longer than typical mt labels.
    if len(hg) >= 4 and hg[0].isalpha() and hg[0].isupper():
        tail = hg[1:]
        if any(ch.isdigit() for ch in tail) and all(ch.isalnum() for ch in hg):
            return "y"

    mt_prefixes = tuple("LMNRHVJTUKIWXABCDFGYZPQSOab")
    if hg[0] in mt_prefixes:
        return "mt"

    return ""


def _finalize_system(df: pd.DataFrame, sheet_name: Optional[Union[str, int]]) -> pd.DataFrame:
    df = df.copy()
    sheet_system = _infer_system_from_sheet_name(sheet_name)

    if sheet_system:
        mask = df["system"].isna() | (df["system"].astype(str).str.strip() == "")
        df.loc[mask, "system"] = sheet_system
        # For dedicated Y / mt sheets, trust the sheet label for all rows.
        df["system"] = df["system"].replace("", sheet_system).fillna(sheet_system)

    def _infer_row(row) -> str:
        existing = str(row.get("system", "") or "").strip().lower()
        if existing in ("mt", "y"):
            return existing
        return _infer_system_from_hg(str(row.get("haplogroup", "") or ""))

    df["system"] = df.apply(_infer_row, axis=1)
    return df


def _clean_rows(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df = df[df["haplogroup"].notna()]
    df = df[df["haplogroup"].astype(str).str.strip() != ""]
    df = df[df["vip_name"].astype(str).str.strip() != ""]
    return df.reset_index(drop=True)


def load_vip_table(
    file_path: str | Path,
    sheet_name: Optional[Union[str, int]] = 0,
) -> pd.DataFrame:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"VIP file not found: {file_path}")

    if file_path.suffix.lower() not in [".xlsx", ".xls"]:
        df = pd.read_csv(file_path)
        df = _normalize_vip_columns(df)
        df = _finalize_system(df, sheet_name=None)
        return _clean_rows(df)

    if sheet_name in (None, "all", "ALL", "*"):
        sheets = pd.read_excel(file_path, sheet_name=None)
        parts = []
        for sname, sdf in sheets.items():
            sdf = _normalize_vip_columns(sdf)
            sdf = _finalize_system(sdf, sheet_name=sname)
            sdf["sheet_name"] = sname
            parts.append(_clean_rows(sdf))
        if not parts:
            return pd.DataFrame(columns=["vip_name", "haplogroup", "system", "source", "note", "category", "sheet_name"])
        return pd.concat(parts, ignore_index=True)

    df = pd.read_excel(file_path, sheet_name=sheet_name)
    df = _normalize_vip_columns(df)
    df = _finalize_system(df, sheet_name=sheet_name)
    return _clean_rows(df)
