# -*- coding: utf-8 -*-

"""The graphical part of a Thermochemistry step"""

import pprint  # noqa: F401
import tkinter as tk
import tkinter.ttk as ttk

from .thermochemistry_parameters import ThermochemistryParameters
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_widgets as sw


class TkThermochemistry(seamm.TkNode):
    """
    The graphical part of a Thermochemistry step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in Thermochemistry_parameters.py

    See Also
    --------
    Thermochemistry, TkThermochemistry,
    ThermochemistryParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        namespace="org.molssi.seamm.tk",
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.namespace = namespace
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )
        self.create_dialog()

    def create_dialog(self):
        """
        Create the dialog. A set of widgets will be chosen by default
        based on what is specified in the Thermochemistry_parameters
        module.

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkThermochemistry.reset_dialog
        """

        frame = super().create_dialog(title="Thermochemistry", widget="notebook")
        # make it large!
        screen_w = self.dialog.winfo_screenwidth()
        screen_h = self.dialog.winfo_screenheight()
        w = int(0.9 * screen_w)
        h = int(0.8 * screen_h)
        x = int(0.05 * screen_w / 2)
        y = int(0.1 * screen_h / 2)

        self.dialog.geometry(f"{w}x{h}+{x}+{y}")

        # Add a frame for the flowchart
        notebook = self["notebook"]
        flowchart_frame = ttk.Frame(notebook)
        self["flowchart frame"] = flowchart_frame
        notebook.add(flowchart_frame, text="Flowchart", sticky=tk.NSEW)

        self.tk_subflowchart = seamm.TkFlowchart(
            master=flowchart_frame,
            flowchart=self.node.subflowchart,
            namespace=self.namespace,
        )
        self.tk_subflowchart.draw()

        # Fill in the control parameters
        # Shortcut for parameters
        P = self.node.parameters

        # thermochemistry frame to isolate widgets
        frame = self["thermochemistry frame"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Thermochemistry Parameters",
            labelanchor="n",
            padding=10,
        )

        for key in ThermochemistryParameters.parameters:
            if key not in ("results",):
                self[key] = P[key].widget(frame)

        # Shown when the Model Chemistry is chosen but none precedes this step
        self["model chemistry note"] = ttk.Label(
            frame,
            text=(
                "Add a Model Chemistry step before this step: it defines the "
                "method (e.g. MOPAC PM6-ORG, or an ORCA DFT model) used for the "
                "energy and the Hessian. Or choose the subflowchart."
            ),
            foreground="red",
            wraplength=500,
            justify=tk.LEFT,
        )

        # Comboboxes whose value changes the layout re-lay out the dialog
        for key in ("engine", "approach"):
            widget = self[key]
            combobox = widget.combobox if hasattr(widget, "combobox") else widget
            combobox.bind("<<ComboboxSelected>>", self.reset_dialog)
            combobox.bind("<Return>", self.reset_dialog)
            combobox.bind("<FocusOut>", self.reset_dialog)

        # and lay them out
        self.reset_dialog()

        self.setup_results()

    def reset_dialog(self, widget=None):
        """Layout the widgets in the dialog.

        The widgets are chosen by default from the information in
        Thermochemistry parameters.

        This function simply lays them out row by row with
        aligned labels. You may wish a more complicated layout that
        is controlled by values of some of the control parameters.
        If so, edit or override this method

        Parameters
        ----------
        widget : Tk Widget = None

        Returns
        -------
        None

        See Also
        --------
        TkThermochemistry.create_dialog
        """

        # Remove any widgets previously packed
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0

        self["thermochemistry frame"].grid(row=row, column=0, sticky=tk.EW, pady=10)
        row += 1
        self.reset_thermochemistry_frame()

        frame.columnconfigure(0, weight=1)

        return row

    def reset_thermochemistry_frame(self, widget=None):
        """Layout the widgets in the thermochemistry frame
        as needed for the current state"""

        approach = self["approach"].get()

        frame = self["thermochemistry frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        # Main controls
        row = 0
        widgets = []
        for key in ("engine", "approach"):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.W)
            widgets.append(self[key])
            row += 1

        # The subflowchart tab applies only to the subflowchart engine
        engine = self["engine"].get()
        notebook = self["notebook"]
        flowchart_frame = self["flowchart frame"]
        if engine == "subflowchart":
            notebook.add(flowchart_frame, text="Flowchart", sticky=tk.NSEW)
        else:
            notebook.hide(flowchart_frame)
            if not self._upstream_has_model_chemistry():
                self["model chemistry note"].grid(
                    row=row, column=0, columnspan=2, sticky=tk.W, pady=(0, 6)
                )
                row += 1

        if approach == "Harmonic approximation":
            for key in (
                "# imaginary modes",
                "spin multiplicity",
                "symmetry number",
                "step size",
                "T",
                "P",
                "print frequencies",
            ):
                self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                widgets.append(self[key])
                row += 1

            if engine == "subflowchart":
                for key in (
                    "on success",
                    "on error",
                ):
                    self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                    widgets.append(self[key])
                    row += 1

        sw.align_labels(widgets, sticky=tk.E)
        # w1 = sw.align_labels(widgets, sticky=tk.E)
        # w2 = sw.align_labels(widgets2, sticky=tk.E)
        # frame.columnconfigure(0, minsize=w1 - w2 + 50)
        frame.columnconfigure(1, weight=1)

    def _upstream_has_model_chemistry(self):
        """True if a Model Chemistry step precedes this one in the flowchart.
        On any error (e.g. the node is not yet linked in) returns False, so the
        reminder is shown."""
        try:
            return any(
                type(node).__name__ == "ModelChemistry"
                and type(node).__module__.startswith("model_chemistry_step")
                for node in self.previous_nodes()
            )
        except Exception:
            return False

    def edit(self):
        """Re-lay out the dialog when it is opened, so the reminder tracks
        the flowchart."""
        if self.dialog is not None:
            self.reset_dialog()
        super().edit()

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkThermochemistry.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
