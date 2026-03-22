"""
analysis_vip.py
---------------
VIP (notable individuals) matching analysis.

This module mirrors the structure of :mod:`analysis_origin` but focuses
entirely on linking the query haplogroup to historical VIP individuals.
It is run alongside the origin analyser; their outputs are merged by
:func:`service.combine_results`.

For Y-DNA queries, VIP haplogroups are resolved to tree nodes before matching
so that relationships are determined via tree traversal rather than string
prefix heuristics alone.

Classes
~~~~~~~
HaplogroupAnalysisResult  — dataclass (identical to the one in analysis_origin).
HaplogroupAnalyzer        — stateful analyser; call ``analyze()`` to run.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Set

import pandas as pd

from vip_matcher import match_vips_for_target
from y_mapper import (
    build_y_mapping_from_aadr,
    convert_tree_labels_to_aadr_terminal,
    load_y_aliases,
    normalize_y_label,
    resolve_y_label_for_tree,
)


@dataclass
class HaplogroupAnalysisResult:
    target: str
    system: str
    resolved_target_for_tree: Optional[str]
    upstream_label: Optional[str]
    downstream_labels: List[str]
    included_labels: List[str]
    matched_samples: pd.DataFrame
    oldest_sample: pd.DataFrame
    early_samples: pd.DataFrame
    country_summary: pd.DataFrame
    candidate_origin_country: Optional[str]
    vip_matches: pd.DataFrame


class HaplogroupAnalyzer:
    def __init__(
        self,
        aadr_df: pd.DataFrame,
        mt_tree,
        y_tree,
        vip_df: Optional[pd.DataFrame] = None,
        y_alias_path: str = "../data/y_aliases.tsv",
    ) -> None:
        self.aadr_df = aadr_df.copy()
        self.mt_tree = mt_tree
        self.y_tree = y_tree
        self.vip_df = vip_df
        self.y_alias_df = load_y_aliases(y_alias_path)

    def get_one_level_upstream_label(self, target: str, tree) -> Optional[str]:
        if target and tree.has_node(target):
            return tree.get_parent(target)
        return None

    def analyze(self, target: str, system: str, early_n: int = 5) -> HaplogroupAnalysisResult:
        system = str(system).strip().lower()
        target = str(target).strip()

        if system not in {"mt", "y"}:
            raise ValueError("system must be 'mt' or 'y'")

        if system == "mt":
            return self._analyze_mt(target=target, early_n=early_n)
        return self._analyze_y(target=target, early_n=early_n)

    # ---------------- mtDNA ----------------
    def _analyze_mt(self, target: str, early_n: int = 5) -> HaplogroupAnalysisResult:
        tree = self.mt_tree
        hg_column = "mt_haplogroup"

        resolved_target_for_tree = target
        upstream_label = self.get_one_level_upstream_label(target, tree)
        downstream_labels = tree.get_descendants(target) if tree.has_node(target) else []

        included_labels = {target}
        if upstream_label:
            included_labels.add(upstream_label)
        included_labels.update(downstream_labels)

        matched_samples = self.filter_aadr_samples(
            df=self.aadr_df,
            haplogroup_column=hg_column,
            included_labels=included_labels,
        )
        matched_samples = self.sort_samples_by_date(matched_samples)
        oldest_sample = self.get_oldest_sample(matched_samples)
        early_samples = self.get_early_samples(matched_samples, n=early_n)
        country_summary = self.summarize_early_countries(early_samples)
        candidate_origin_country = self.infer_candidate_origin_country(country_summary)

        vip_matches = self.match_vips_mt(target=target, tree=tree)

        return HaplogroupAnalysisResult(
            target=target,
            system="mt",
            resolved_target_for_tree=resolved_target_for_tree,
            upstream_label=upstream_label,
            downstream_labels=downstream_labels,
            included_labels=sorted(included_labels),
            matched_samples=matched_samples,
            oldest_sample=oldest_sample,
            early_samples=early_samples,
            country_summary=country_summary,
            candidate_origin_country=candidate_origin_country,
            vip_matches=vip_matches,
        )

    def match_vips_mt(self, target: str, tree) -> pd.DataFrame:
        if self.vip_df is None:
            return self._empty_vip_result()

        return match_vips_for_target(
            vip_df=self.vip_df,
            target=target,
            system="mt",
            tree=tree,
            include_exact=True,
            include_upstream=True,
            include_downstream=True,
            match_column="haplogroup",
            output_column="haplogroup",
            keep_related_without_direct_relation=True,
        )

    # ---------------- Y chromosome ----------------
    def _analyze_y(self, target: str, early_n: int = 5) -> HaplogroupAnalysisResult:
        tree = self.y_tree
        terminal_to_isogg, isogg_to_terminal = build_y_mapping_from_aadr(self.aadr_df)

        resolved_target_for_tree = resolve_y_label_for_tree(
            label=target,
            tree=tree,
            terminal_to_isogg=terminal_to_isogg,
            alias_df=self.y_alias_df,
        )

        upstream_label: Optional[str] = None
        downstream_labels: List[str] = []
        included_labels: Set[str] = set()

        if resolved_target_for_tree is None:
            included_labels.add(normalize_y_label(target))
        else:
            upstream_label = self.get_one_level_upstream_label(resolved_target_for_tree, tree)
            downstream_labels = tree.get_descendants(resolved_target_for_tree) if tree.has_node(resolved_target_for_tree) else []

            labels_in_tree = {resolved_target_for_tree}
            if upstream_label:
                labels_in_tree.add(upstream_label)
            labels_in_tree.update(downstream_labels)

            included_labels = convert_tree_labels_to_aadr_terminal(
                labels_in_tree=labels_in_tree,
                isogg_to_terminal=isogg_to_terminal,
                alias_df=self.y_alias_df,
            )
            included_labels.add(normalize_y_label(target))

        matched_samples = self.filter_aadr_y_samples(self.aadr_df, included_labels)
        matched_samples = self.sort_samples_by_date(matched_samples)
        oldest_sample = self.get_oldest_sample(matched_samples)
        early_samples = self.get_early_samples(matched_samples, n=early_n)
        country_summary = self.summarize_early_countries(early_samples)
        candidate_origin_country = self.infer_candidate_origin_country(country_summary)

        vip_matches = self.match_vips_y(
            resolved_target_for_tree=resolved_target_for_tree,
            tree=tree,
            terminal_to_isogg=terminal_to_isogg,
        )

        return HaplogroupAnalysisResult(
            target=target,
            system="y",
            resolved_target_for_tree=resolved_target_for_tree,
            upstream_label=upstream_label,
            downstream_labels=downstream_labels,
            included_labels=sorted(included_labels),
            matched_samples=matched_samples,
            oldest_sample=oldest_sample,
            early_samples=early_samples,
            country_summary=country_summary,
            candidate_origin_country=candidate_origin_country,
            vip_matches=vip_matches,
        )

    def match_vips_y(
        self,
        resolved_target_for_tree: Optional[str],
        tree,
        terminal_to_isogg: dict,
    ) -> pd.DataFrame:
        if self.vip_df is None or resolved_target_for_tree is None:
            return self._empty_vip_result()

        vip_df = self.vip_df.copy()
        vip_df["original_haplogroup"] = vip_df["haplogroup"].astype(str).str.strip().apply(normalize_y_label)
        vip_df["resolved_haplogroup"] = vip_df["haplogroup"].apply(
            lambda x: resolve_y_label_for_tree(
                label=x,
                tree=tree,
                terminal_to_isogg=terminal_to_isogg,
                alias_df=self.y_alias_df,
            )
        )

        vip_df = vip_df[vip_df["original_haplogroup"] != ""].copy()

        # For Y, keep only downstream VIPs.
        # Use the original ISOGG-style haplogroup for relationship checks so
        # deep downstream clades are not collapsed into the same resolved ancestor.
        return match_vips_for_target(
            vip_df=vip_df,
            target=normalize_y_label(resolved_target_for_tree),
            system="y",
            tree=tree,
            include_exact=False,
            include_upstream=False,
            include_downstream=True,
            match_column="original_haplogroup",
            output_column="original_haplogroup",
            keep_related_without_direct_relation=False,
        )

    # ---------------- shared ----------------
    def _empty_vip_result(self) -> pd.DataFrame:
        return pd.DataFrame(
            columns=[
                "target",
                "system",
                "vip_name",
                "vip_haplogroup",
                "vip_resolved_haplogroup",
                "relation",
                "common_ancestor",
                "source",
                "note",
            ]
        )

    def filter_aadr_samples(
        self,
        df: pd.DataFrame,
        haplogroup_column: str,
        included_labels: Set[str],
    ) -> pd.DataFrame:
        out = df.copy()
        out[haplogroup_column] = out[haplogroup_column].astype(str).str.strip()
        out = out[out[haplogroup_column].isin(included_labels)]
        return out.reset_index(drop=True)

    def filter_aadr_y_samples(self, df: pd.DataFrame, included_labels: Set[str]) -> pd.DataFrame:
        norm_df = df.copy()
        norm_df["_y_terminal"] = norm_df["y_haplogroup"].apply(normalize_y_label)
        norm_df["_y_isogg"] = norm_df["y_haplogroup_isogg"].apply(normalize_y_label)
        labels = {normalize_y_label(x) for x in included_labels if normalize_y_label(x)}
        mask = norm_df["_y_terminal"].isin(labels) | norm_df["_y_isogg"].isin(labels)
        return norm_df.loc[mask].drop(columns=["_y_terminal", "_y_isogg"]).reset_index(drop=True)

    def sort_samples_by_date(self, df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df.copy()

        out = df.copy()
        if "date_mean_bp" in out.columns:
            out["date_mean_bp"] = pd.to_numeric(out["date_mean_bp"], errors="coerce")
            out = out.sort_values(by="date_mean_bp", ascending=False, na_position="last")
        return out.reset_index(drop=True)

    def get_oldest_sample(self, df: pd.DataFrame) -> pd.DataFrame:
        if df.empty:
            return df.copy()
        return df.iloc[[0]].copy()

    def get_early_samples(self, df: pd.DataFrame, n: int = 5) -> pd.DataFrame:
        if df.empty:
            return df.copy()
        return df.head(n).copy()

    def summarize_early_countries(self, df: pd.DataFrame) -> pd.DataFrame:
        if df.empty or "political_entity" not in df.columns:
            return pd.DataFrame(columns=["political_entity", "sample_count", "oldest_bp"])

        out = df.copy()
        out["political_entity"] = out["political_entity"].fillna("Unknown")
        if "date_mean_bp" in out.columns:
            out["date_mean_bp"] = pd.to_numeric(out["date_mean_bp"], errors="coerce")

        summary = (
            out.groupby("political_entity", dropna=False)
            .agg(
                sample_count=("political_entity", "size"),
                oldest_bp=("date_mean_bp", "max"),
            )
            .reset_index()
            .sort_values(["sample_count", "oldest_bp", "political_entity"], ascending=[False, False, True])
            .reset_index(drop=True)
        )
        return summary

    def infer_candidate_origin_country(self, country_summary: pd.DataFrame) -> Optional[str]:
        if country_summary.empty:
            return None
        return str(country_summary.iloc[0]["political_entity"])
