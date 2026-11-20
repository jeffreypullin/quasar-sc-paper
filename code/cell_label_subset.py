
from __future__ import annotations

import pandas as pd

# When several datasets are run together, this dict is the union of combined-type
# definitions across all of them. Constraint: a combined-type name must map to a
# single global member set (i.e. names must be globally unique across datasets).
# onek1k defines the entries below; tremor currently adds none. If a future
# dataset needs to redefine one of these names, pass the dataset through the
# matching scripts and key this dict by (dataset, cell_type) instead.
COMBINED_MEMBERS: dict[str, tuple[str, ...]] = {
    "B_all": ("B IN", "B Mem"),
    "T_all": ("CD4 NC", "CD4 ET", "CD4 SOX4", "CD8 ET", "CD8 NC", "CD8 S100B"),
    # Cytoimmgen: cell_label is "<Compartment> <timepoint>" after harmonise.
    # Update membership if timepoints in params.dataset_configs.cytoimmgen change.
    "Naive_all": ("Naive 0h", "Naive 16h", "Naive 40h", "Naive 64h"),
    "Memory_all": ("Memory 0h", "Memory 16h", "Memory 40h", "Memory 64h"),
}


def cell_label_mask(cell_label_series: pd.Series, cell_type: str) -> pd.Series:
    if cell_type in COMBINED_MEMBERS:
        return cell_label_series.isin(COMBINED_MEMBERS[cell_type])
    return cell_label_series == cell_type.replace("_", " ")
