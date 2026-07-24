"""Plain string-patching helpers for a py3Dmol `_make_html()` scene, vendored
verbatim from `dd_viewer.scene` (only these four functions -- the rest of
that module builds a single-receptor/pose scene via rdkit, which dd_idea
doesn't need; see `dd_idea/viewer3d/__init__.py`).
"""
from __future__ import annotations

import json
import re
from typing import Optional

_VIEWER_RENDER_CALL = re.compile(r"viewer_(\d+)\.render\(\);")


def get_viewer_variable(html: str) -> Optional[str]:
    """The py3Dmol-generated JS variable name (e.g. `"viewer_123456"`) that
    a `_make_html()` string's scene is loaded into, or `None` if `html`
    doesn't look like a py3Dmol scene at all. Each `_make_html()` call
    picks a fresh random suffix, so this has to be re-extracted per HTML
    string rather than assumed constant -- used by both
    `html_with_camera_events`/`html_with_initial_view` to know which
    variable to patch calls onto.
    """
    match = _VIEWER_RENDER_CALL.search(html)
    return f"viewer_{match.group(1)}" if match else None


def html_with_initial_view(html: str, view: list) -> str:
    """Patch a py3Dmol `_make_html()` string to apply a saved camera view (a
    `getView()`-shaped array) right after the scene renders, instead of
    leaving the scene's own default `zoomTo()` auto-fit as the final word.
    A no-op if `html` doesn't look like a py3Dmol scene, or if `view` is
    falsy (`None`/`[]` -- no saved view to restore).
    """
    if not view:
        return html
    viewer_var = get_viewer_variable(html)
    if viewer_var is None:
        return html
    match = _VIEWER_RENDER_CALL.search(html)
    snippet = f"\n{viewer_var}.setView({json.dumps(view)});\n"
    return html.replace(match.group(0), match.group(0) + snippet, 1)


def html_with_camera_events(html: str) -> str:
    """Patch a py3Dmol `_make_html()` string to report its camera state and
    render-readiness to the parent window via postMessage, for the
    `dd_idea_3d` double-buffered component (`dd_idea.viewer3d.component`):

    - On every drag/zoom/touch, posts `{plviewerCameraUpdate: true, view:
      [...]}` (the 3Dmol `getView()` array). The component's own JS keeps
      the latest one in memory (its execution context survives Streamlit
      reruns, unlike this scene's) and re-applies it to each newly-loaded
      scene -- which is what makes the camera position survive a widget
      interaction instead of snapping back to the default zoomTo fit every
      time. The *same* report is also sent once, proactively, right after
      the initial paint (see below) -- without this, a scene the user never
      manually dragged/zoomed leaves the component with nothing saved to
      restore, so the next widget interaction (e.g. a checkbox that changes
      which structures are shown) would fall through to that new scene's
      own default `zoomTo()` fit instead, visibly reframing the camera even
      though the user never touched it.
    - Two animation frames after the initial render (once the paint has
      actually landed), posts `{plviewerReady: true}` -- the signal the
      component waits for before swapping this scene into view, so updates
      read as a cross-fade instead of a flash to blank.

    Both are harmless no-ops if nothing is listening (e.g. a plain
    `st.iframe` embed, or Jupyter).
    """
    viewer_var = get_viewer_variable(html)
    if viewer_var is None:
        return html
    match = _VIEWER_RENDER_CALL.search(html)
    snippet = f"""
try {{
  var __plvSave = function() {{
    try {{ parent.postMessage({{plviewerCameraUpdate: true, view: {viewer_var}.getView()}}, "*"); }} catch (e) {{}}
  }};
  document.addEventListener('mouseup', __plvSave);
  document.addEventListener('wheel', __plvSave, {{passive: true}});
  document.addEventListener('touchend', __plvSave);
}} catch (e) {{}}
requestAnimationFrame(function() {{
  requestAnimationFrame(function() {{
    try {{ __plvSave(); }} catch (e) {{}}
    try {{ parent.postMessage({{plviewerReady: true}}, "*"); }} catch (e) {{}}
  }});
}});
"""
    return html.replace(match.group(0), match.group(0) + snippet, 1)


_FIXED_SIZE_DIV_STYLE = re.compile(r'style="position: relative; width: \d+px; height: \d+px;"')


def html_fill_container(html: str) -> str:
    """Patch a py3Dmol `_make_html()` string so its viewer div fills 100% of
    whatever element embeds it, instead of the fixed pixel `width`/`height`
    baked in from the `width=`/`height=` passed to `py3Dmol.view(...)`.
    """
    html = _FIXED_SIZE_DIV_STYLE.sub('style="position: relative; width: 100%; height: 100%;"', html, count=1)
    reset = "<style>html, body { margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; }</style>"
    return reset + html
