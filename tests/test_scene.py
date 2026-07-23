"""Fast, offline unit tests for scene.py's pure-Python helpers. Not
covered: build_overlay_view's actual rendered py3Dmol/WebGL output, per
this project family's manual-only convention for that side of things
(see dd_seqalign's own test_scene.py)."""
from dd_compare.scene import PALETTE, PDB_PALETTE, _readable_font_color


def test_readable_font_color_picks_black_on_light_backgrounds():
    assert _readable_font_color("#ffe119") == "black"  # PDB_PALETTE's bright yellow
    assert _readable_font_color("#bcf60c") == "black"  # PDB_PALETTE's bright lime
    assert _readable_font_color("#ffffff") == "black"


def test_readable_font_color_picks_white_on_dark_backgrounds():
    assert _readable_font_color("#444444") == "white"  # REFERENCE_COLOR
    assert _readable_font_color("#1f77b4") == "white"  # PALETTE's blue
    assert _readable_font_color("#000000") == "white"


def test_readable_font_color_covers_every_palette_entry():
    # Every currently-shipped palette color should resolve to one of the
    # two valid choices without raising -- guards against a malformed hex
    # string slipping into either palette unnoticed.
    for color in [*PALETTE, *PDB_PALETTE, "#444444"]:
        assert _readable_font_color(color) in ("black", "white")
