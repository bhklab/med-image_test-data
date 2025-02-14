from collections import defaultdict
import json
import pandas as pd
# snakefile to download data
configfile: "nbia.yaml"
COLLECTION_NAMES = list(config["datasets"].keys())
rule all:
  input:
    expand("results/{collection}.tar.gz", collection=COLLECTION_NAMES),
    expand("results/summaries/{collection}_summary.md", collection=COLLECTION_NAMES),
    expand("metadata/{collection}.json", collection=COLLECTION_NAMES)
  output:
    report = "results/combined_summary.md",
  run:
        
    all_collection_modality_sets = defaultdict(list)
    # get all the modalities for each collection
    for collection in COLLECTION_NAMES:
      with open(f"metadata/{collection}.json") as metadata_file:
        metadata = json.load(metadata_file)
        for item in metadata:
          all_collection_modality_sets[collection].append(item["Modality"])
    # we want to have a dataframe where columns are UNIQUE modalities
    # and rows are collections
    # we will use this to generate a markdown table

    # get all unique modalities
    all_modalities = set()
    for modalities in all_collection_modality_sets.values():
      all_modalities.update(modalities)
    
    # sort the modalitiese
    # CT, MR, SEG, RTSTRUCT, RTDOSE, RTPLAN come first, then everything else
    sorted_modalities = sorted(all_modalities, key=lambda x: (x not in ["CT", "MR", "SEG", "RTSTRUCT", "RTDOSE", "RTPLAN"], x))
    cols = ["CT", "MR", "SEG", "RTSTRUCT", "RTDOSE", "RTPLAN"] + [x for x in sorted_modalities if x not in ["CT", "MR", "SEG", "RTSTRUCT", "RTDOSE", "RTPLAN"]]

    # create a dataframe with all modalities as columns
    # and collections as rows
    df = pd.DataFrame(index=COLLECTION_NAMES, columns=cols)

    for collection, modalities in all_collection_modality_sets.items():
      for modality in modalities:
        df.at[collection, modality] = all_collection_modality_sets[collection].count(modality)

    # sort rows by collection name
    df.sort_index(inplace=True)
    # replace NaN with empty string
    df.fillna("", inplace=True)

    print(df)


    # use python to concatenate all the summary files, then sort by collection name
    summaries = []
    for collection in COLLECTION_NAMES:
      with open(f"results/{collection}_summary.md") as summary:
        summaries.append((collection, summary.read()))
    
    summaries.sort(key=lambda x: x[0])
    
    with open("results/combined_summary.md", "a") as combined_summary:
      # write the modalities table
      combined_summary.write("## Modalities Summary\n\n")
      combined_summary.write(df.to_markdown())

    
      combined_summary.write("\n\n")
      for collection, content in summaries:
        combined_summary.write(content)
        combined_summary.write("\n")



rule compress:
  input:
    collection_dir = "procdata/{collection}"
  output:
    gzip_dir = "results/{collection}.tar.gz"
  shell:
    # this command will first tar the directory
    # uses -C to change to the directory
    # uses -cf to create a new tar file
    # uses --exclude to exclude the .snakemake_timestamp file (because we use the 'directory' function in snakemake)
    # then pipes the output to pigz to compress it
    """
    tar --exclude='.snakemake_timestamp' -C procdata -cf - {wildcards.collection} | pigz -9 -c > {output.gzip_dir}
    """

rule summarize_metadata:
  input:
    metadata_file = "metadata/{collection}.json"
  output:
    summary_file = "results/summaries/{collection}_summary.md"
  script:
    "workflow/scripts/summarize_metadata.py"

rule download_collection:
  output:
    collection_dir = directory("procdata/{collection}"),
    metadata_file = "metadata/{collection}.json"
  retries: 3
  log:
    "logs/download_collection/download_{collection}.log"
  script:
    "workflow/scripts/download_collection.py"