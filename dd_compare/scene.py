"""py3Dmol multi-protein overlay scene, built from `structalign`'s
superposed AlphaFold-model coordinate files. Adapted from
`dd_seqalign.scene.build_overlay_view` for the cross-protein case: that
function highlights the *same* canonical positions across multiple
structures of one protein; here every protein has its own numbering, so
each structure gets its own list of residues to highlight (translated
position-by-position by `pocketmap.py`). No ligand handling is needed here
(every structure is an AlphaFold model, always apo), which is the other
simplification relative to `dd_seqalign.scene`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Sequence, Union

import py3Dmol

PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]
REFERENCE_COLOR = "#444444"
SITE_COLOR = "yellow"

# 3Dmol.js's built-in "*Carbon" colorschemes (e.g. "yellowCarbon") tint carbon
# and leave every other element at its RasMol default, but only for a fixed
# set of named CSS colors -- this is the same RasMol table, copied here so
# any per-protein hex can be used as the carbon tint via a
# `{"prop": "elem", "map": {...}}` colorscheme instead of a scheme name.
_HETERO_ELEMENT_COLORS = {
    "H": "#ffffff", "He": "#ffc0cb", "Li": "#b22222", "B": "#00ff00",
    "N": "#8f8fff", "O": "#f00000", "F": "#daa520", "Na": "#0000ff",
    "Mg": "#228b22", "Al": "#808090", "Si": "#daa520", "P": "#ffa500",
    "S": "#ffc832", "Cl": "#00ff00", "Ca": "#808090",
}


def _carbon_tint_scheme(carbon_color: str) -> dict:
    return {"prop": "elem", "map": {**_HETERO_ELEMENT_COLORS, "C": carbon_color}}


def assign_colors(labels: Sequence[str], reference_label: str) -> Dict[str, str]:
    """One color per label, cycling through `PALETTE`; the reference always
    gets the same neutral gray so it reads as "the reference", not just
    another protein in the cycle."""
    colors: Dict[str, str] = {}
    i = 0
    for label in labels:
        if label == reference_label:
            colors[label] = REFERENCE_COLOR
        else:
            colors[label] = PALETTE[i % len(PALETTE)]
            i += 1
    return colors


def build_overlay_view(
    structures: Sequence[dict], reference_label: str, *,
    colors: Optional[Dict[str, str]] = None, width: Union[int, str] = "100%", height: Union[int, str] = 600,
) -> py3Dmol.view:
    """`structures`: each a dict with `label` (UniProt accession),
    `pdb_path` (superposed coordinates, from `structalign.align_structures`'s
    `aligned_pdb`), `chain_id`, and `site_resseqs` (that protein's *own*
    AFDB residue numbers to highlight -- see `pocketmap.py`; empty/omitted
    draws no highlight for that structure, e.g. positions with no
    counterpart in this particular protein)."""
    view = py3Dmol.view(width=width, height=height)
    colors = colors or assign_colors([s["label"] for s in structures], reference_label)

    for model_index, s in enumerate(structures):
        pdb_text = Path(s["pdb_path"]).read_text()
        view.addModel(pdb_text, "pdb")
        color = colors.get(s["label"], "gray")
        chain_sel = {"model": model_index, "chain": s["chain_id"]}

        # Every chain is explicitly styled to nothing first, since 3Dmol.js
        # falls back to a default line/wireframe rendering for any atom
        # left unstyled rather than hiding it.
        view.setStyle({"model": model_index}, {})
        view.setStyle(chain_sel, {"cartoon": {"color": color}})

        site = s.get("site_resseqs")
        if site:
            view.addStyle(
                {"model": model_index, "chain": s["chain_id"], "resi": list(site)},
                {"stick": {"colorscheme": _carbon_tint_scheme(SITE_COLOR), "radius": 0.25}},
            )

    view.zoomTo()
    return view
