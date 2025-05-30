#!/usr/bin/env python3

from collections import defaultdict
from pathlib import Path
import pandas as pd
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import snakemake  # type: ignore

METADATA_FILES = [Path(f) for f in snakemake.input.metadata_files]
SUMMARY_FILES = [Path(f) for f in snakemake.input.summaries]
output_path = Path(snakemake.output.combined_summary)
combined_metadata = Path(snakemake.output.combined_metadata)
all_collection_modality_sets: dict[str, list[str]] = defaultdict(list)

# Re-discover collection list (glob resets on full iteration)
for csv_file in METADATA_FILES:
    collection = csv_file.stem
    df = pd.read_csv(csv_file)
    if "Modality" not in df.columns:
        raise ValueError(f"{csv_file} missing 'Modality' column")
    all_collection_modality_sets[collection] = df["Modality"].tolist()

# Get all unique modalities
all_modalities = set()
for modalities in all_collection_modality_sets.values():
    all_modalities.update(modalities)

priority_modalities = ["CT", "MR", "SEG", "RTSTRUCT", "RTDOSE", "RTPLAN"]
sorted_modalities = sorted(
    all_modalities, key=lambda x: (x not in priority_modalities, x)
)
cols = priority_modalities + [
    m for m in sorted_modalities if m not in priority_modalities
]

# Create DataFrame
modality_counts = defaultdict(int)
seen = set()

df = pd.DataFrame(index=all_collection_modality_sets.keys(), columns=cols)
# drop duplicates
df = df.drop_duplicates()

for collection, modalities in all_collection_modality_sets.items():
    for modality in modalities:
        count = modalities.count(modality)
        df.at[collection, modality] = count

        if (collection, modality) not in seen:
            modality_counts[modality] += count
            seen.add((collection, modality))

df.sort_index(inplace=True)
df.fillna("", inplace=True)

for modality in modality_counts:
    if modality in df.columns:
        df.at["Total", modality] = modality_counts[modality]
    else:
        df.at["Total", modality] = 0


# Append summaries
summaries = []
for collection in all_collection_modality_sets:
    summary_path = Path(f"results/summaries/{collection}_summary.md")
    if summary_path.exists():
        summaries.append((collection, summary_path.read_text()))

summaries.sort(key=lambda x: x[0])


with output_path.open("w") as f:
    f.write("## Modalities Summary\n\n")
    f.write(df.to_markdown())
    f.write("\n\n")
    for collection, content in summaries:
        f.write(content + "\n")

# combine all csvs into one
df = pd.concat([pd.read_csv(f) for f in METADATA_FILES], ignore_index=True, sort=False)
df = df.drop_duplicates()
df.to_csv(combined_metadata, index=False)
