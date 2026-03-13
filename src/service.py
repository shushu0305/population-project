from types import SimpleNamespace

from aadr_parser import load_aadr_table
from analysis_origin import HaplogroupAnalyzer as OriginHaplogroupAnalyzer
from analysis_vip import HaplogroupAnalyzer as VipHaplogroupAnalyzer
from tree_parser import PhyloTree
from vip_parser import load_vip_table


def load_trees(mt_tree_path: str, y_tree_path: str):
    mt_tree = PhyloTree(mt_tree_path)
    mt_tree.load()

    y_tree = PhyloTree(y_tree_path)
    y_tree.load()

    return mt_tree, y_tree


def combine_results(origin_result, vip_result):
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


def run_haplogroup_analysis(
    system: str,
    target: str,
    aadr_path: str,
    vip_path: str | None = None,
    vip_sheet: str | int | None = None,
    mt_tree_path: str = "../data/trees/mt_phyloTree_b17_Tree2.txt",
    y_tree_path: str = "../data/trees/chrY_hGrpTree_isogg2016.txt",
    early_n: int = 5,
):
    mt_tree, y_tree = load_trees(mt_tree_path, y_tree_path)

    aadr_df = load_aadr_table(aadr_path, sep="\t")

    vip_df = None
    if vip_path:
        if isinstance(vip_sheet, str) and vip_sheet.isdigit():
            vip_sheet = int(vip_sheet)

        if vip_sheet is None:
            vip_sheet = "Y" if system == "y" else "mtDNA"

        vip_df = load_vip_table(vip_path, sheet_name=vip_sheet)

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
        target=target,
        system=system,
        early_n=early_n,
    )
    vip_result = vip_analyzer.analyze(
        target=target,
        system=system,
        early_n=early_n,
    )

    return combine_results(origin_result, vip_result)
