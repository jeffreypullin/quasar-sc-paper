
from __future__ import annotations

import pandas as pd

COMBINED_MEMBERS: dict[str, tuple[str, ...]] = {
    "B_all": ("B IN", "B Mem"),
    "CD4_T_all": ("CD4 NC", "CD4 ET", "CD4 SOX4"),
    "CD8_T_all": ("CD8 ET", "CD8 NC", "CD8 S100B"),
}


def cell_label_mask(cell_label_series: pd.Series, cell_type: str) -> pd.Series:
    if cell_type in COMBINED_MEMBERS:
        return cell_label_series.isin(COMBINED_MEMBERS[cell_type])
    return cell_label_series == cell_type
