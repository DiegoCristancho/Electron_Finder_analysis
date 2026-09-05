"""Run the preserved legacy calculation with an injected chain and ROOT module."""

import ast
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
LEGACY = PROJECT / "OLD" / "other" / "Distances_Between_cluster_track.py"


def assigns(node, name):
    return isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == name for target in node.targets)


def run_reference(events, root):
    """Execute the original helpers, histograms, event loop, and derived counts.

    Only original file loading and plotting/printing are excluded. No formulas
    or selections are translated here, so this provides an independent oracle.
    """
    tree = ast.parse(LEGACY.read_text())
    chain_start = next(i for i, node in enumerate(tree.body) if assigns(node, "chain"))
    histogram_start = next(i for i, node in enumerate(tree.body) if assigns(node, "h_d3D_noEoP"))
    report_start = next(i for i, node in enumerate(tree.body) if i > histogram_start
                        and isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                        and isinstance(node.value.func, ast.Name) and node.value.func.id == "print")
    prefix = [node for node in tree.body[:chain_start]
              if not (isinstance(node, ast.Import) and any(alias.name == "ROOT" for alias in node.names))]
    core = ast.Module(body=prefix + tree.body[histogram_start:report_start], type_ignores=[])
    namespace = {"ROOT": root, "chain": events}
    exec(compile(core, str(LEGACY), "exec"), namespace)
    return namespace
