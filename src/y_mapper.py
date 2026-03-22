"""
y_mapper.py
-----------
Y-DNA label resolution utilities.

Y-chromosome haplogroup labels exist in two incompatible formats:

* **Terminal SNP format** (used in AADR):  e.g. ``R-M269``, ``I-M253``
* **ISOGG hierarchical format** (used in the 2016 tree):  e.g. ``R1b1a1``, ``I1``

This module bridges the two formats through a four-step resolution pipeline
(see :func:`resolve_y_label_for_tree`) and provides helpers for normalisation
and bidirectional mapping.

Key functions
~~~~~~~~~~~~~
normalize_y_label(label)                        — strip noise, return clean str
load_y_aliases(path)                            — load the alias TSV table
build_y_mapping_from_aadr(aadr_df)             — build terminal↔ISOGG dicts
resolve_y_label_for_tree(label, tree, ...)      — resolve any Y label to tree node
convert_tree_labels_to_aadr_terminal(labels, ...)— reverse: tree → AADR format
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd


def normalize_y_label(label: str) -> str:
    if label is None:
        return ""

    label = str(label).strip()
    label = label.replace("?", "").replace("*", "")
    label = label.replace(" ", "")
    if "(" in label:
        label = label.split("(")[0].strip()
    return label


def load_y_aliases(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=["raw_label", "aadr_isogg", "tree_label"])

    df = pd.read_csv(path, sep="\t")

    for col in ["raw_label", "aadr_isogg", "tree_label"]:
        if col not in df.columns:
            df[col] = None
        df[col] = df[col].apply(normalize_y_label)

    return df


def build_y_mapping_from_aadr(aadr_df: pd.DataFrame) -> tuple[dict, dict]:
    required_cols = ["y_haplogroup", "y_haplogroup_isogg"]
    for col in required_cols:
        if col not in aadr_df.columns:
            raise ValueError(f"Missing required AADR column: {col}")

    df = aadr_df[required_cols].copy()
    df["y_haplogroup"] = df["y_haplogroup"].apply(normalize_y_label)
    df["y_haplogroup_isogg"] = df["y_haplogroup_isogg"].apply(normalize_y_label)

    df = df[
        (df["y_haplogroup"] != "") &
        (df["y_haplogroup_isogg"] != "")
    ].drop_duplicates()

    terminal_to_isogg = {}
    isogg_to_terminal = {}

    for _, row in df.iterrows():
        terminal = row["y_haplogroup"]
        isogg = row["y_haplogroup_isogg"]
        terminal_to_isogg.setdefault(terminal, isogg)
        isogg_to_terminal.setdefault(isogg, terminal)

    return terminal_to_isogg, isogg_to_terminal


def _resolve_from_aliases(label: str, tree, alias_df: Optional[pd.DataFrame]) -> Optional[str]:
    if alias_df is None or alias_df.empty:
        return None

    hit = alias_df[alias_df["raw_label"] == label]
    if not hit.empty:
        tree_label = normalize_y_label(hit.iloc[0]["tree_label"])
        if tree.has_node(tree_label):
            return tree_label

    hit = alias_df[alias_df["aadr_isogg"] == label]
    if not hit.empty:
        tree_label = normalize_y_label(hit.iloc[0]["tree_label"])
        if tree.has_node(tree_label):
            return tree_label

    return None


def _nearest_existing_prefix(label: str, tree) -> Optional[str]:
    """
    Fallback for ISOGG-style labels such as R1b1a1b1a1a2b1.
    If the exact node is absent in the tree, walk upward by trimming the last
    character until an existing ancestor node is found.
    """
    label = normalize_y_label(label)
    if not label:
        return None

    if tree.has_node(label):
        return label

    candidate = label
    while len(candidate) > 1:
        candidate = candidate[:-1]
        if tree.has_node(candidate):
            return candidate

    return None


def resolve_y_label_for_tree(
    label: str,
    tree,
    terminal_to_isogg: dict,
    alias_df: Optional[pd.DataFrame] = None,
) -> Optional[str]:
    """
    Resolve a Y label to a tree-usable label.

    Priority:
    1. label itself in tree
    2. alias raw_label / aadr_isogg -> tree_label
    3. AADR terminal -> isogg, if isogg in tree
    4. nearest existing prefix ancestor in tree
    """
    label = normalize_y_label(label)

    if not label:
        return None

    if tree.has_node(label):
        return label

    alias_hit = _resolve_from_aliases(label, tree=tree, alias_df=alias_df)
    if alias_hit:
        return alias_hit

    mapped_isogg = terminal_to_isogg.get(label)
    if mapped_isogg:
        mapped_isogg = normalize_y_label(mapped_isogg)
        if tree.has_node(mapped_isogg):
            return mapped_isogg

        alias_hit = _resolve_from_aliases(mapped_isogg, tree=tree, alias_df=alias_df)
        if alias_hit:
            return alias_hit

        prefix_hit = _nearest_existing_prefix(mapped_isogg, tree)
        if prefix_hit:
            return prefix_hit

    prefix_hit = _nearest_existing_prefix(label, tree)
    if prefix_hit:
        return prefix_hit

    return None


def convert_tree_labels_to_aadr_terminal(
    labels_in_tree: set[str],
    isogg_to_terminal: dict,
    alias_df: Optional[pd.DataFrame] = None,
) -> set[str]:
    out = set()

    alias_tree_to_raw = {}
    if alias_df is not None and not alias_df.empty:
        for _, row in alias_df.iterrows():
            raw_label = normalize_y_label(row["raw_label"])
            tree_label = normalize_y_label(row["tree_label"])
            if raw_label and tree_label:
                alias_tree_to_raw.setdefault(tree_label, raw_label)

    for label in labels_in_tree:
        label = normalize_y_label(label)
        if label in alias_tree_to_raw:
            out.add(alias_tree_to_raw[label])
        elif label in isogg_to_terminal:
            out.add(isogg_to_terminal[label])
        else:
            out.add(label)

    return out
