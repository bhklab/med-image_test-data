from collections import defaultdict
import json
import pandas as pd

import os

ACCESS_TYPE = os.environ.get("ACCESS_TYPE", "public").upper()
print(f"ACCESS_TYPE: {ACCESS_TYPE}")

if ACCESS_TYPE == "PUBLIC":
  configfile: "nbia_datasets.yaml"
  configfile: "nbia_datasets.yaml"
elif ACCESS_TYPE == "PRIVATE":
  configfile: "nbia_datasets_private.yaml"

  NBIA_USERNAME = os.environ.get("NBIA_USERNAME", None)
  NBIA_PASSWORD =  os.environ.get("NBIA_PASSWORD", None)
  if NBIA_USERNAME is None or NBIA_PASSWORD is None:
    raise ValueError("NBIA_USERNAME and NBIA_PASSWORD must be set for private access.")

  os.environ["LOGIN__NBIA_USERNAME"] = NBIA_USERNAME
  os.environ["LOGIN__NBIA_PASSWORD"] = NBIA_PASSWORD

elif ACCESS_TYPE == "BOTH":
  import yaml
  # read both yaml files 
  with open("nbia_datasets.yaml") as public_file:
    public_config = yaml.safe_load(public_file)

  with open("nbia_datasets_private.yaml") as private_file:
    private_config = yaml.safe_load(private_file)

  # merge the two configs
  config = public_config
  config["datasets"].update(private_config["datasets"])
else:
  raise ValueError(f"Invalid ACCESS_TYPE: {ACCESS_TYPE}. Must be 'public' or 'private'.")

COLLECTION_NAMES = list(config["datasets"].keys())
print(f"COLLECTION_NAMES: {COLLECTION_NAMES}")


rule all:
  input:
    collection_data = expand("results/{collection}.tar.gz", collection=COLLECTION_NAMES),
    summaries = expand("results/summaries/{collection}_summary.md", collection=COLLECTION_NAMES),
    metadata_files = expand("metadata/{collection}.csv", collection=COLLECTION_NAMES)
  output:
    combined_summary="results/combined_summary.md",
    combined_metadata=f"results/{ACCESS_TYPE}-combined_metadata.csv",
  script:
    "workflow/scripts/create_combined_summary.py"

rule compress:
  input:
    collection_dir = "procdata/{collection}"
  output:
    gzip_dir = "results/{collection}.tar.gz"
  shell:
    """
    INPUT_DIR=$(dirname {input.collection_dir})
    tar --exclude='.snakemake_timestamp' -C $INPUT_DIR -cf - {wildcards.collection} | pigz -9 -c > {output.gzip_dir}
    """
    # this command will first tar the directory
    # uses -C to change to the directory
    # uses -cf to create a new tar file
    # uses --exclude to exclude the .snakemake_timestamp file (because we use the 'directory' function in snakemake)
    # then pipes the output to pigz to compress it

rule summarize_metadata:
  input:
    metadata_file = "metadata/{collection}.csv"
  output:
    summary_file = "results/summaries/{collection}_summary.md"
  script:
    "workflow/scripts/summarize_metadata.py"

rule download_collection:
  input:
    metadata_file = "metadata/{collection}.csv"
  input:
    metadata_file = "metadata/{collection}.csv"
  output:
    collection_dir = directory("procdata/{collection}"),
  log:
    "logs/download_collection/download_{collection}.log"
  retries: 3
  script:
    "workflow/scripts/download_collection.py"

rule get_metadata:
  output:
    metadata_file = "metadata/{collection}.csv"
  log:
    "logs/get_metadata/{collection}.log"
  script:
    "workflow/scripts/get_metadata.py"


# rule all:
#   input:
#     expand("results/{collection}.tar.gz", collection=COLLECTION_NAMES),
#     expand("results/summaries/{collection}_summary.md", collection=COLLECTION_NAMES),
#     expand("metadata/{collection}.json", collection=COLLECTION_NAMES)
#   output:
#     report = "results/combined_summary.md",
#   run:
        
#     all_collection_modality_sets = defaultdict(list)
#     # get all the modalities for each collection
#     for collection in COLLECTION_NAMES:
#       with open(f"metadata/{collection}.json") as metadata_file:
#         metadata = json.load(metadata_file)
#         for item in metadata:
#           all_collection_modality_sets[collection].append(item["Modality"])
#     # we want to have a dataframe where columns are UNIQUE modalities
#     # and rows are collections
#     # we will use this to generate a markdown table

#     # get all unique modalities
#     all_modalities = set()
#     for modalities in all_collection_modality_sets.values():
#       all_modalities.update(modalities)
    
#     # sort the modalitiese
#     # CT, MR, SEG, RTSTRUCT, RTDOSE, RTPLAN come first, then everything else
#     sorted_modalities = sorted(all_modalities, key=lambda x: (x not in ["CT", "MR", "SEG", "RTSTRUCT", "RTDOSE", "RTPLAN"], x))
#     cols = ["CT", "MR", "SEG", "RTSTRUCT", "RTDOSE", "RTPLAN"] + [x for x in sorted_modalities if x not in ["CT", "MR", "SEG", "RTSTRUCT", "RTDOSE", "RTPLAN"]]

#     # create a dataframe with all modalities as columns
#     # and collections as rows
#     df = pd.DataFrame(index=COLLECTION_NAMES, columns=cols)

#     for collection, modalities in all_collection_modality_sets.items():
#       for modality in modalities:
#         df.at[collection, modality] = all_collection_modality_sets[collection].count(modality)

#     # sort rows by collection name
#     df.sort_index(inplace=True)
#     # replace NaN with empty string
#     df.fillna("", inplace=True)

#     # use python to concatenate all the summary files, then sort by collection name
#     summaries = []
#     for collection in COLLECTION_NAMES:
#       with open(f"results/summaries/{collection}_summary.md") as summary:
#         summaries.append((collection, summary.read()))
    
#     summaries.sort(key=lambda x: x[0])
    
#     with open("results/combined_summary.md", "a") as combined_summary:
#       # write the modalities table
#       combined_summary.write("## Modalities Summary\n\n")
#       combined_summary.write(df.to_markdown())

    
#       combined_summary.write("\n\n")
#       for collection, content in summaries:
#         combined_summary.write(content)
#         combined_summary.write("\n")

# rule compress:
#   input:
#     collection_dir = "procdata/sorted/{collection}"
#   output:
#     gzip_dir = "results/{collection}.tar.gz"
#   shell:
#     """
#     INPUT_DIR=$(dirname {input.collection_dir})
#     tar --exclude='.snakemake_timestamp' -C $INPUT_DIR -cf - {wildcards.collection} | pigz -9 -c > {output.gzip_dir}
#     """
#     # this command will first tar the directory
#     # uses -C to change to the directory
#     # uses -cf to create a new tar file
#     # uses --exclude to exclude the .snakemake_timestamp file (because we use the 'directory' function in snakemake)
#     # then pipes the output to pigz to compress it

# rule dicomsort:
#   input: 
#     collection_dir = "procdata/unsorted/{collection}"
#   output:
#     sorted_dir = directory("procdata/sorted/{collection}")
#   shell:
#     """
#     imgtools dicomsort \
#     --action move \
#     {input.collection_dir} \
#     {output.sorted_dir}/%PatientID/%Modality_Series-%SeriesInstanceUID/

#     rm -rf {input.collection_dir}
#     """

# rule summarize_metadata:
#   input:
#     metadata_file = "metadata/{collection}.json"
#   output:
#     summary_file = "results/summaries/{collection}_summary.md"
#   script:
#     "workflow/scripts/summarize_metadata.py"

# rule download_collection:
#   output:
#     collection_dir = directory("procdata/unsorted/{collection}"),
#     metadata_file = "metadata/{collection}.json"
#   retries: 3
#   log:
#     "logs/download_collection/download_{collection}.log"
#   script:
#     "workflow/scripts/download_collection.py"