from pathlib import Path

import pandas as pd

DATA_PATH = Path("data/raw/GSE226105_rLog-Normalized_counts.csv")

df = pd.read_csv(DATA_PATH)