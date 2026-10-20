"""Isolates excluded for contamination: `contamination_flag == "contaminated"` in screening_map_results.csv
(CheckM contamination above 10 %, set by immensekansasii/scripts/06_screening_map_results.py). One rule for every script."""
from pathlib import Path

import pandas as pd


def contaminated_nr(root: Path) -> set[int]:
    m = pd.read_csv(root / "data/imm/screening_map_results.csv", usecols=["NR", "contamination_flag"])
    return {int(n) for n in m.loc[m.contamination_flag == "contaminated", "NR"]}
