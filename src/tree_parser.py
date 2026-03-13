from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, List, Optional, Set


@dataclass(frozen=True)
class Edge:
    child: str
    parent: str


class PhyloTree:
    """
    Parse a simple haplogroup tree file in `child parent` format.

    Supported examples from your data:
    - chrY_hGrpTree_isogg2016.txt
    - mt_phyloTree_b17_Tree2.txt
    """

    def __init__(self, file_path: str | Path, root_names: Optional[Set[str]] = None) -> None:
        self.file_path = Path(file_path)
        self.root_names = root_names or {"Root", "mt-MRCA"}

        self.parent_map: Dict[str, str] = {}
        self.children_map: DefaultDict[str, List[str]] = defaultdict(list)
        self.nodes: Set[str] = set()
        self.edges: List[Edge] = []

    def load(self) -> None:
        if not self.file_path.exists():
            raise FileNotFoundError(f"Tree file not found: {self.file_path}")

        with self.file_path.open("r", encoding="utf-8") as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue

                parts = line.split()
                if len(parts) < 2:
                    continue

                child, parent = parts[0], parts[1]

                if self._is_header_or_invalid(child, parent):
                    continue

                self.parent_map[child] = parent
                self.children_map[parent].append(child)
                self.nodes.add(child)
                self.nodes.add(parent)
                self.edges.append(Edge(child=child, parent=parent))

    def _is_header_or_invalid(self, child: str, parent: str) -> bool:
        header_pairs = {
            ("Root", "#"),
            ("mt-MRCA", "#"),
        }
        if (child, parent) in header_pairs:
            return True

        # Skip obvious column headers if present
        lowered = {child.lower(), parent.lower()}
        if "root" in lowered and "#" in {child, parent}:
            return True

        return False

    def has_node(self, node: str) -> bool:
        return node in self.nodes

    def get_parent(self, node: str) -> Optional[str]:
        return self.parent_map.get(node)

    def get_children(self, node: str) -> List[str]:
        return list(self.children_map.get(node, []))

    def get_ancestors(self, node: str, include_self: bool = False) -> List[str]:
        if not self.has_node(node):
            return []

        ancestors: List[str] = [node] if include_self else []
        current = node

        while current in self.parent_map:
            parent = self.parent_map[current]
            ancestors.append(parent)
            current = parent

            if parent in self.root_names:
                break

        return ancestors

    def get_descendants(self, node: str, include_self: bool = False) -> List[str]:
        if not self.has_node(node):
            return []

        descendants: List[str] = [node] if include_self else []
        stack = list(self.get_children(node))

        while stack:
            current = stack.pop()
            descendants.append(current)
            stack.extend(self.get_children(current))

        return descendants

    def get_lineage_to_root(self, node: str) -> List[str]:
        """
        Return lineage as:
        [node, parent, grandparent, ..., root]
        """
        return self.get_ancestors(node, include_self=True)

    def get_depth(self, node: str) -> int:
        return max(0, len(self.get_ancestors(node)) - 1)

    def is_ancestor_of(self, ancestor: str, node: str) -> bool:
        return ancestor in self.get_ancestors(node)

    def is_descendant_of(self, node: str, ancestor: str) -> bool:
        return ancestor in self.get_ancestors(node)

    def relationship(self, node_a: str, node_b: str) -> str:
        """
        Return one of:
        - exact
        - ancestor
        - descendant
        - unrelated
        - missing
        """
        if not self.has_node(node_a) or not self.has_node(node_b):
            return "missing"
        if node_a == node_b:
            return "exact"
        if self.is_ancestor_of(node_a, node_b):
            return "ancestor"
        if self.is_descendant_of(node_a, node_b):
            return "descendant"
        return "unrelated"

    def print_subtree(self, node: str, max_depth: Optional[int] = None) -> None:
        if not self.has_node(node):
            print(f"[missing] {node}")
            return

        def _print(current: str, level: int) -> None:
            print("  " * level + current)
            if max_depth is not None and level >= max_depth:
                return
            for child in sorted(self.get_children(current)):
                _print(child, level + 1)

        _print(node, 0)

    def bfs_subtree(self, node: str) -> List[str]:
        if not self.has_node(node):
            return []

        result: List[str] = []
        queue: deque[str] = deque([node])

        while queue:
            current = queue.popleft()
            result.append(current)
            for child in self.get_children(current):
                queue.append(child)

        return result

