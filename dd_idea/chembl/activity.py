"""Counting ChEMBL bioactivity data for a target -- how much SAR data
exists, not the data itself (see `dd_chembl` for actually pulling and
modeling it). Cheap by design: ChEMBL's `/activity.json` reports
`page_meta.total_count` even with `limit=1`, so getting a count never
means paginating through every record.
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from typing import Tuple

CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_ACTIVITY_SEARCH = CHEMBL_API_BASE + "/activity.json"

# Same default as dd_chembl/dd_chembl/fetch.py: binding assays only, with a
# pChEMBL value (i.e. actually quantitatively comparable, not every raw
# assay result ChEMBL has on file for the target).
DEFAULT_ASSAY_TYPES: Tuple[str, ...] = ("B",)


def count_activities(target_chembl_id: str, *, assay_types: Tuple[str, ...] = DEFAULT_ASSAY_TYPES) -> int:
    """How many ChEMBL activity records exist for `target_chembl_id`,
    restricted to `assay_types` (default: binding, 'B') and requiring a
    pChEMBL value -- the same filter `dd_chembl.fetch.fetch_activities`
    uses, so this count matches what dd_chembl would actually pull for
    QSAR training."""
    query = urllib.parse.urlencode({
        "target_chembl_id": target_chembl_id,
        "assay_type__in": ",".join(assay_types),
        "pchembl_value__isnull": "false",
        "limit": 1,
        "format": "json",
    })
    with urllib.request.urlopen(f"{CHEMBL_ACTIVITY_SEARCH}?{query}") as fh:
        body = json.load(fh)
    return body["page_meta"]["total_count"]
