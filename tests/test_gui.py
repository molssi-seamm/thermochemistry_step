# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialog: create it and re-lay it out for every choice
that drives the layout. Skipped when no display is available."""

import pytest


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


def _mapped(widget):
    return widget.grid_info() != {}


def test_dialog_layouts(tk_node):
    tk_node.create_dialog()
    tk_node.reset_dialog()
    notebook = tk_node["notebook"]
    flowchart_frame = tk_node["flowchart frame"]

    # The default is the Model Chemistry: no subflowchart tab, no subdirectory
    # handling, and a reminder since no Model Chemistry precedes the step
    assert tk_node["engine"].get() == "Model Chemistry"
    assert notebook.tab(flowchart_frame, "state") == "hidden"
    assert not _mapped(tk_node["on success"])
    assert _mapped(tk_node["model chemistry note"])
    assert _mapped(tk_node["symmetry number"])

    # The subflowchart: the tab is back, with the subdirectory handling
    tk_node["engine"].set("subflowchart")
    tk_node.reset_dialog()
    assert notebook.tab(flowchart_frame, "state") == "normal"
    assert _mapped(tk_node["on success"])
    assert not _mapped(tk_node["model chemistry note"])

    tk_node["engine"].set("Model Chemistry")
    tk_node.reset_dialog()
    assert notebook.tab(flowchart_frame, "state") == "hidden"

    for approach in tk_node["approach"].combobox.cget("values"):
        tk_node["approach"].set(approach)
        tk_node.reset_dialog()


def test_cbh_layout(tk_node):
    tk_node.create_dialog()
    tk_node["approach"].set("Isodesmic (CBH) with ATcT references")
    tk_node.reset_dialog()
    # Computed energies, harmonic: the optimization and step size, no T or P,
    # no property names
    assert _mapped(tk_node["cbh rungs"]) and _mapped(tk_node["ATcT version"])
    assert _mapped(tk_node["fragment convergence"])
    assert _mapped(tk_node["step size"])
    assert not _mapped(tk_node["T"]) and not _mapped(tk_node["P"])
    assert not _mapped(tk_node["energy property"])
    assert not _mapped(tk_node["on success"])

    assert _mapped(tk_node["fragment store file"])
    tk_node["fragment store"].set("ignore")
    tk_node.reset_dialog()
    assert not _mapped(tk_node["fragment store file"])

    # No thermal corrections (MOPAC): no step size
    tk_node["thermal corrections"].set("none (the energy is the enthalpy, e.g. MOPAC)")
    tk_node.reset_dialog()
    assert not _mapped(tk_node["step size"])
    assert _mapped(tk_node["fragment convergence"])

    # No optimization: no convergence settings
    tk_node["optimize fragments"].set("no")
    tk_node.reset_dialog()
    assert not _mapped(tk_node["fragment convergence"])
    assert not _mapped(tk_node["fragment max steps"])

    # Existing results: the property names, the ZPE/enthalpy ones only when
    # harmonic
    tk_node["cbh energies"].set("existing results (run nothing)")
    tk_node.reset_dialog()
    assert _mapped(tk_node["energy property"])
    assert not _mapped(tk_node["ZPE property"])
    assert not _mapped(tk_node["optimize fragments"])
    tk_node["thermal corrections"].set("harmonic")
    tk_node.reset_dialog()
    assert _mapped(tk_node["ZPE property"])
    assert _mapped(tk_node["enthalpy property"])
    assert not _mapped(tk_node["step size"])

    # The subflowchart engine adds its error handling for computed energies
    tk_node["cbh energies"].set("finite-difference harmonic")
    tk_node["engine"].set("subflowchart")
    tk_node.reset_dialog()
    assert _mapped(tk_node["on success"])
    assert tk_node["notebook"].tab(tk_node["flowchart frame"], "state") == "normal"

    # Back to the harmonic approach: nothing of the CBH left
    tk_node["approach"].set("Harmonic approximation")
    tk_node.reset_dialog()
    assert not _mapped(tk_node["cbh rungs"])
    assert _mapped(tk_node["T"])
