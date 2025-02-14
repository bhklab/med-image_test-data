import json
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import snakemake  # type: ignore

# Load JSON data
json_file = Path(snakemake.input["metadata_file"])
output_file = Path(snakemake.output["summary_file"])
collection_name = snakemake.wildcards.collection
with open(json_file, "r") as f:
    data = json.load(f)

# Aggregate data
summary = defaultdict(lambda: defaultdict(int))
for entry in data:
    patient_id = entry["PatientID"]
    modality = entry["Modality"]
    summary[patient_id][modality] += 1

# Generate Markdown table
header = (
    f"# {collection_name} Summary\n\n"
    "| Patient ID           | Modality  | Number of Series    |\n"
    "|----------------------|-----------|---------------------|\n"
)

rows = []

for patient_id, modalities in summary.items():
    # write the patient_id in the first row
    # the patient id in the other rows will be empty
    for i, (modality, count) in enumerate(modalities.items()):
        if i == 0:
            rows.append(f"| {patient_id:<20} | {modality:<9} | {count:<19} |")
        else:
            rows.append(f"| {'':<20} | {modality:<9} | {count:<19} |")



# Combine header and rows
markdown_output = header + "\n".join(rows) + "\n"

# Save Markdown output
output_file.write_text(markdown_output)
print(f"Summary written to {output_file}")
