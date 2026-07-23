"""py3Dmol multi-protein overlay scene, built from `structalign`'s
superposed AlphaFold-model coordinate files, plus zero or more independent
real-RCSB-structure overlays per protein (see `pdbstruct.py`/`pipeline.py`'s
`pdb_structures` report field). Adapted from
`dd_seqalign.scene.build_overlay_view` for the cross-protein case: that
function highlights the *same* canonical positions across multiple
structures of one protein; here every protein has its own numbering, so
each structure gets its own list of residues to highlight (translated
position-by-position by `pocketmap.py`).
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
PDB_CARTOON_OPACITY = 0.6  # visible-but-distinct from the solid AlphaFold cartoon (0.35 was too faint to read)
LIGAND_COLOR = "#ff9900"  # fallback when no per-ligand color was assigned

# Distinct from PALETTE (protein cartoon colors) so a ligand's stick color
# never visually collides with any protein's own cartoon, and multiple
# ligands shown together (even across different proteins) stay
# distinguishable from each other.
LIGAND_PALETTE = [
    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
]

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
    another protein in the cycle. Duplicate labels (e.g. a protein's
    AlphaFold entry and its own real-PDB-overlay entries, which
    deliberately share a label so they render in the same color) collapse
    to one color, not one each."""
    colors: Dict[str, str] = {}
    i = 0
    for label in labels:
        if label in colors:
            continue
        if label == reference_label:
            colors[label] = REFERENCE_COLOR
        else:
            colors[label] = PALETTE[i % len(PALETTE)]
            i += 1
    return colors


def assign_ligand_colors(keys: Sequence[str]) -> Dict[str, str]:
    """One color per key (e.g. an `f"{accession}:{pdb_id}"` string),
    cycling through `LIGAND_PALETTE`. Assign this over the *full* set of
    available real-structure entries (not just the ones currently checked
    to show) so a given ligand's color stays stable as the user toggles
    other structures on/off."""
    return {key: LIGAND_PALETTE[i % len(LIGAND_PALETTE)] for i, key in enumerate(keys)}


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
    """`structures`: a flat list of independent entries to draw -- normally
    one `kind="afdb"` entry per shown protein (its AlphaFold model) plus
    zero or more `kind="pdb"` entries (real RCSB structures, already
    superposed onto the reference frame -- see `pdbstruct.py`/
    `pipeline.py`'s `pdb_structures` report field). Keeping every real
    structure as its own list entry, rather than nesting them under a
    parent protein entry, lets the caller show/hide a protein's AlphaFold
    model and each of its real structures independently of one another.

    Each entry: `label` (used to look up a shared per-protein `color` from
    `colors` when `color` isn't given directly on the entry -- an afdb
    entry and its protein's own pdb entries normally share a label so they
    render in the same color), `pdb_path`, `chain_id`, `site_resseqs`/
    `site_labels` (that structure's *own* numbering -- for `kind="pdb"` see
    `pipeline.py`'s `pocket_resseq` mapping, since a real PDB entry's
    residue numbers don't equal canonical UniProt position, unlike an
    AlphaFold model; `site_labels` is an optional `{resnum: "K33"}` map,
    used only when `label_residues` is True). `kind` (`"afdb"`, the
    default: solid cartoon; `"pdb"`: semi-transparent cartoon, see
    `PDB_CARTOON_OPACITY`, plus thinner pocket-residue sticks). For
    `kind="pdb"` entries only: `ligand_resname` (drawn as sticks if given
    -- the actual payoff of overlaying a real structure at all) and
    `ligand_color` (that ligand's own distinct color, e.g. from
    `assign_ligand_colors`, so multiple ligands shown together stay
    visually distinguishable from each other; falls back to the fixed
    `LIGAND_COLOR` if not given)."""
    view = py3Dmol.view(width=width, height=height)
    colors = colors or assign_colors([s["label"] for s in structures], reference_label)

    for model_index, s in enumerate(structures):
        pdb_text = Path(s["pdb_path"]).read_text()
        view.addModel(pdb_text, "pdb")
        color = s.get("color") or colors.get(s["label"], "gray")
        kind = s.get("kind", "afdb")
        opacity = PDB_CARTOON_OPACITY if kind == "pdb" else 1.0
        chain_sel = {"model": model_index, "chain": s["chain_id"]}

        # Every chain is explicitly styled to nothing first, since 3Dmol.js
        # falls back to a default line/wireframe rendering for any atom
        # left unstyled rather than hiding it.
        view.setStyle({"model": model_index}, {})
        view.setStyle(chain_sel, {"cartoon": {"color": color, "opacity": opacity}})

        site = s.get("site_resseqs")
        if site:
            stick_radius = 0.18 if kind == "pdb" else 0.25
            view.addStyle(
                {"model": model_index, "chain": s["chain_id"], "resi": list(site)},
                {"stick": {"colorscheme": _carbon_tint_scheme(SITE_COLOR), "radius": stick_radius}},
            )
            if label_residues:
                site_labels = {r: (s.get("site_labels") or {}).get(r, str(r)) for r in site}
                _add_site_labels(view, s["pdb_path"], s["chain_id"], site_labels, color)

        if kind == "pdb":
            ligand_resname = s.get("ligand_resname")
            if ligand_resname:
                ligand_color = s.get("ligand_color") or LIGAND_COLOR
                view.addStyle(
                    {"model": model_index, "resn": ligand_resname},
                    {"stick": {"colorscheme": _carbon_tint_scheme(ligand_color), "radius": 0.3}},
                )

    view.zoomTo()
    return view
