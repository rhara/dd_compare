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
from typing import Dict, Optional, Sequence, Tuple, Union

import py3Dmol

PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]
REFERENCE_COLOR = "#444444"
SITE_COLOR = "yellow"
PDB_CARTOON_OPACITY = 0.35  # thinner/more transparent than the AlphaFold cartoon, so the two stay visually distinct
LIGAND_COLOR = "#ff9900"

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


def _ca_coords(pdb_path: str, chain_id: str, resnums: Sequence[int]) -> Dict[int, Tuple[float, float, float]]:
    """CA atom position for each of `resnums` in the given chain, read
    directly from the PDB file text (fixed-column parsing, same convention
    as `pdbio.py`) -- py3Dmol/3Dmol.js has no synchronous "read back a
    loaded model's coordinates" API from Python, so label placement needs
    its own tiny parse of the same file `addModel` was given."""
    wanted = set(resnums)
    coords: Dict[int, Tuple[float, float, float]] = {}
    for line in Path(pdb_path).read_text().splitlines():
        if line[:6] != "ATOM  " or line[12:16].strip() != "CA" or line[21] != chain_id:
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        if resseq in wanted and resseq not in coords:
            coords[resseq] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return coords


def _add_site_labels(view: py3Dmol.view, pdb_path: str, chain_id: str, site_labels: Dict[int, str], color: str) -> None:
    for resnum, coord in _ca_coords(pdb_path, chain_id, list(site_labels)).items():
        view.addLabel(site_labels.get(resnum, str(resnum)), {
            "position": {"x": coord[0], "y": coord[1], "z": coord[2]},
            "backgroundColor": color, "backgroundOpacity": 0.8,
            "fontColor": "white", "fontSize": 11, "showBackground": True, "borderThickness": 0.4,
        })


def build_overlay_view(
    structures: Sequence[dict], reference_label: str, *,
    colors: Optional[Dict[str, str]] = None, width: Union[int, str] = "100%", height: Union[int, str] = 600,
    label_residues: bool = False,
) -> py3Dmol.view:
    """`structures`: each a dict with `label` (UniProt accession),
    `pdb_path` (superposed coordinates, from `structalign.align_structures`'s
    `aligned_pdb`), `chain_id`, `site_resseqs` (that protein's *own*
    AFDB residue numbers to highlight -- see `pocketmap.py`; empty/omitted
    draws no highlight for that structure, e.g. positions with no
    counterpart in this particular protein), `site_labels` (optional
    `{resnum: "K33"}` text for each of `site_resseqs`, used only when
    `label_residues` is True -- one 3D text marker per lining residue,
    anchored at its CA position and colored like that protein's cartoon),
    and an optional `pdb_overlay` dict (from `pdbstruct.py`/`pipeline.py`'s
    `pdb` report field) -- a real RCSB structure for the same protein,
    already superposed onto the reference frame, drawn as a second,
    semi-transparent cartoon alongside the AlphaFold one (never replacing
    it -- see README "Why always AlphaFold"): `pdb_path`, `chain_id`,
    `site_resseqs`/`site_labels` in *that structure's own* numbering (see
    `pipeline.py`'s `pocket_resseq` mapping -- a real PDB entry's residue
    numbers don't equal canonical UniProt position, unlike an AlphaFold
    model), and `ligand_resname` (the bound ligand to render as sticks, if
    any -- the actual payoff of overlaying a real structure at all)."""
    view = py3Dmol.view(width=width, height=height)
    colors = colors or assign_colors([s["label"] for s in structures], reference_label)

    next_model = 0
    for s in structures:
        model_index = next_model
        next_model += 1
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
            if label_residues:
                site_labels = {r: (s.get("site_labels") or {}).get(r, str(r)) for r in site}
                _add_site_labels(view, s["pdb_path"], s["chain_id"], site_labels, color)

        pdb_overlay = s.get("pdb_overlay")
        if pdb_overlay:
            pdb_model_index = next_model
            next_model += 1
            pdb_overlay_text = Path(pdb_overlay["pdb_path"]).read_text()
            view.addModel(pdb_overlay_text, "pdb")
            pdb_chain_sel = {"model": pdb_model_index, "chain": pdb_overlay["chain_id"]}

            view.setStyle({"model": pdb_model_index}, {})
            view.setStyle(pdb_chain_sel, {"cartoon": {"color": color, "opacity": PDB_CARTOON_OPACITY}})

            pdb_site = pdb_overlay.get("site_resseqs")
            if pdb_site:
                view.addStyle(
                    {"model": pdb_model_index, "chain": pdb_overlay["chain_id"], "resi": list(pdb_site)},
                    {"stick": {"colorscheme": _carbon_tint_scheme(SITE_COLOR), "radius": 0.18}},
                )
                if label_residues:
                    pdb_site_labels = {r: (pdb_overlay.get("site_labels") or {}).get(r, str(r)) for r in pdb_site}
                    _add_site_labels(view, pdb_overlay["pdb_path"], pdb_overlay["chain_id"], pdb_site_labels, color)

            ligand_resname = pdb_overlay.get("ligand_resname")
            if ligand_resname:
                view.addStyle(
                    {"model": pdb_model_index, "resn": ligand_resname},
                    {"stick": {"colorscheme": _carbon_tint_scheme(LIGAND_COLOR), "radius": 0.3}},
                )

    view.zoomTo()
    return view
