"""
config.py
---------
Centralised path configuration for Haplogroup Discover.

All data-file paths and runtime defaults are defined here so that the rest
of the codebase never hard-codes paths.  Every value can be overridden at
runtime by setting the corresponding environment variable — no code changes
required.

Environment variables
~~~~~~~~~~~~~~~~~~~~~
AADR_PATH     Path to the AADR annotation file   (default: data/aadr/AADR Annotations 2025.xlsx)
VIP_PATH      Path to the VIP Excel file          (default: data/vip/VIPHaplogroups.xlsx)
MT_TREE_PATH  Path to the mtDNA tree file         (default: data/trees/mt_phyloTree_b17_Tree2.txt)
Y_TREE_PATH   Path to the Y-chromosome tree file  (default: data/trees/chrY_hGrpTree_isogg2016.txt)
EARLY_N       Number of early samples for origin  (default: 5)
"""
from __future__ import annotations

import os
from pathlib import Path

# The level above src/ is the project root directory
_PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _resolve(env_key: str, default_relative: str) -> str:
    """Return a file path, preferring the environment variable if set.

    Parameters
    ----------
    env_key:
        Name of the environment variable to check first.
    default_relative:
        Path relative to the project root, used when the env var is absent.

    Returns
    -------
    str
        Absolute path as a string.
    """
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
