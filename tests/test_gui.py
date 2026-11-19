# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialog of the CBH approach: create it and re-lay it out
for the choices that drive the conformer settings. Skipped when no display is
available."""

import pytest

pytest.importorskip("seamm_thermochemistry")


@pytest.fixture()
def tk_node():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw
    import seamm

    Pmw.initialise(root)
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.tk"
    )
    node = flowchart.create_node("Thermochemistry")
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get("Thermochemistry")
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    yield tk_node
    root.destroy()


def _shown(tk_node, key):
    return tk_node[key].grid_info() != {}


def test_conformers_optimized_only_when_it_applies(tk_node):
    tk_node.create_dialog()
    approaches = [
        a
        for a in tk_node.node.parameters["approach"].enumeration
        if a.startswith("Isodesmic")
    ]
    tk_node["approach"].set(approaches[0])
    tk_node["preset"].set("custom")
    tk_node["cbh energies"].set(tk_node.node.parameters["cbh energies"].enumeration[0])
    tk_node["optimize fragments"].set("yes")
    tk_node["initial conformer"].set("lowest of the embeddings (MMFF)")
    tk_node.reset_dialog()
    assert _shown(tk_node, "conformer embeddings")
    assert _shown(tk_node, "conformers optimized")
    # nothing optimized: the force field's lowest, so no candidates to optimize
    tk_node["optimize fragments"].set("no")
    tk_node.reset_dialog()
    assert not _shown(tk_node, "conformers optimized")
    tk_node["optimize fragments"].set("yes")
    # neither the fragments nor the target from the embeddings
    tk_node["initial conformer"].set("single embedding")
    tk_node["target conformer"].set("as given")
    tk_node.reset_dialog()
    assert not _shown(tk_node, "conformer embeddings")
    assert not _shown(tk_node, "conformers optimized")
    tk_node["target conformer"].set("lowest of the embeddings (MMFF)")
    tk_node.reset_dialog()
    assert _shown(tk_node, "conformers optimized")
