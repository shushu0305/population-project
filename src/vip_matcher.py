from __future__ import annotations

from typing import Optional

import pandas as pd


def normalize_haplogroup_label(label: str) -> str:
    if label is None:
        return ""

    label = str(label).strip()
    for suffix in ["?", "*"]:
        label = label.replace(suffix, "")
    label = label.replace(" ", "")
    if "(" in label:
        label = label.split("(")[0].strip()

    return label


def infer_system_from_label(label: str) -> Optional[str]:
    if not label:
        return None

    label = normalize_haplogroup_label(label)

    if "-" in label:
        return "y"

    mt_prefixes = (
        "L", "M", "N", "R", "H", "V", "J", "T", "U", "K", "I", "W", "X",
        "A", "B", "C", "D", "F", "G", "Y", "Z", "P", "Q", "S", "O"
    )
    if label[:1] in mt_prefixes and "-" not in label:
        return "mt"

    return None


def same_major_clade_mt(target: str, vip_hg: str) -> bool:
    return bool(target and vip_hg and target[0] == vip_hg[0])


def classify_relationship(
    target: str,
    vip_hg: str,
    tree,
    system: str,
) -> str:
    target = normalize_haplogroup_label(target)
    vip_hg = normalize_haplogroup_label(vip_hg)

    if not target or not vip_hg:
        return "missing"

    if target == vip_hg:
        return "exact"

    # For unresolved Y ISOGG labels, fall back to prefix-style lineage checks.
    if system == "y":
        if vip_hg.startswith(target):
            return "downstream"
        if target.startswith(vip_hg):
            return "upstream"

    if not tree.has_node(target) or not tree.has_node(vip_hg):
        return "missing"

    if system == "mt" and not same_major_clade_mt(target, vip_hg):
        return "unrelated"

    ancestors = set(tree.get_ancestors(target))
    descendants = set(tree.get_descendants(target))

    if vip_hg in descendants:
        return "downstream"
    if vip_hg in ancestors:
        return "upstream"

    return "related"


def find_common_ancestor(node_a: str, node_b: str, tree) -> Optional[str]:
    """
    Return the nearest common ancestor of two nodes in the tree.
    For unresolved Y ISOGG-style labels, fall back to longest common prefix.
    """
    node_a = normalize_haplogroup_label(node_a)
    node_b = normalize_haplogroup_label(node_b)

    if tree.has_node(node_a) and tree.has_node(node_b):
        lineage_a = [node_a] + tree.get_ancestors(node_a)
        lineage_b = [node_b] + tree.get_ancestors(node_b)

        lineage_b_set = set(lineage_b)

        for node in lineage_a:
            if node in lineage_b_set:
                return node

    # prefix fallback, mainly for Y labels like R1b1a1b vs R1b1a1b1a1a2b1
    max_len = min(len(node_a), len(node_b))
    prefix = []
    for i in range(max_len):
        if node_a[i] != node_b[i]:
            break
        prefix.append(node_a[i])

    return "".join(prefix) or None


def match_vips_for_target(
    vip_df: pd.DataFrame,
    target: str,
    system: str,
    tree,
    include_upstream: bool = True,
    include_downstream: bool = True,
    include_exact: bool = True,
    match_column: str = "haplogroup",
    output_column: str = "haplogroup",
    keep_related_without_direct_relation: bool = True,
) -> pd.DataFrame:
    """
    For mtDNA:
      returns exact / upstream / downstream / related
    For Y:
      also returns common ancestor if both labels resolve in tree
    """
    target = normalize_haplogroup_label(target)
    system = system.strip().lower()

    results = []

    for _, row in vip_df.iterrows():
        vip_name = row.get("vip_name")
        vip_hg_for_match = normalize_haplogroup_label(row.get(match_column, ""))
        vip_hg_for_output = normalize_haplogroup_label(row.get(output_column, ""))
        vip_system = str(row.get("system", "")).strip().lower()
        source = row.get("source")
        note = row.get("note")

        if not vip_system:
            inferred = infer_system_from_label(vip_hg_for_output)
            vip_system = inferred or ""

        if vip_system != system:
            continue

        relation = classify_relationship(
            target=target,
            vip_hg=vip_hg_for_match,
            tree=tree,
            system=system,
        )

        common_ancestor = find_common_ancestor(target, vip_hg_for_match, tree)

        if relation == "exact" and not include_exact:
            continue
        if relation == "upstream" and not include_upstream:
            continue
        if relation == "downstream" and not include_downstream:
            continue

        if relation in {"missing", "unrelated"}:
            continue

        if relation == "related" and not keep_related_without_direct_relation:
            continue

        results.append(
            {
                "target": target,
                "system": system,
                "vip_name": vip_name,
                "vip_haplogroup": vip_hg_for_output,
                "vip_resolved_haplogroup": vip_hg_for_match,
                "relation": relation,
                "common_ancestor": common_ancestor,
                "source": source,
                "note": note,
            }
        )

    if not results:
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

    out = pd.DataFrame(results)

    relation_order = {
        "exact": 0,
        "downstream": 1,
        "upstream": 2,
        "related": 3,
    }
    out["relation_rank"] = out["relation"].map(relation_order).fillna(99)
    out = out.sort_values(["relation_rank", "vip_name"]).drop(columns=["relation_rank"])

    return out.reset_index(drop=True)