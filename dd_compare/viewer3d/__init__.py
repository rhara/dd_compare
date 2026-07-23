"""Camera-stable py3Dmol embedding for Streamlit, vendored from `dd_viewer`
(itself vendored, unmodified, into every sibling `dd_*` app that needs a 3D
tab) -- trimmed to only the pieces dd_compare actually uses.

The full `dd_viewer` package also ships a single-receptor/single-pose scene
builder (`scene.py`), ligand/interaction detection (`interactions.py`), and
`Receptor`/`Pose` I/O (`io.py`) -- none of which dd_compare needs, and the
first of which imports `rdkit` at module level. dd_compare has no
small-molecule handling of its own (see README "Why no rdkit"), so only the
plain string-patching helpers (`html_with_camera_events`/`html_fill_container`/
`get_viewer_variable`, from `htmlpatch.py` here) and the double-buffered
Streamlit component (`component.py`) are vendored -- importing `dd_compare.viewer3d`
never pulls in rdkit.
"""
from .htmlpatch import get_viewer_variable, html_fill_container, html_with_camera_events, html_with_initial_view

__all__ = [
    "get_viewer_variable",
    "html_fill_container",
    "html_with_camera_events",
    "html_with_initial_view",
    "view3d",
]


def __getattr__(name):
    # Same lazy-import trick as the original `dd_viewer/__init__.py`: importing
    # `dd_compare.viewer3d` at all (e.g. from `pipeline.py`, which never touches
    # Streamlit) shouldn't pull in `streamlit` or trip its "missing
    # ScriptRunContext" warning outside an actual `streamlit run`.
    if name == "view3d":
        from .component import view3d as _view3d
        return _view3d
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
