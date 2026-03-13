from __future__ import annotations

import os
from pathlib import Path

# src/ 的上一级就是项目根目录
_PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _resolve(env_key: str, default_relative: str) -> str:
    """优先读环境变量，否则相对项目根目录拼接绝对路径。"""
    value = os.getenv(env_key)
    if value:
        return value
    return str(_PROJECT_ROOT / default_relative)


DATA_CONFIG: dict = {
    "AADR_PATH":    _resolve("AADR_PATH",    "data/aadr/AADR Annotations 2025.xlsx"),
    "VIP_PATH":     _resolve("VIP_PATH",     "data/vip/VIPHaplogroups.xlsx"),
    "MT_TREE_PATH": _resolve("MT_TREE_PATH", "data/trees/mt_phyloTree_b17_Tree2.txt"),
    "Y_TREE_PATH":  _resolve("Y_TREE_PATH",  "data/trees/chrY_hGrpTree_isogg2016.txt"),
    "EARLY_N":      int(os.getenv("EARLY_N", "5")),
}
