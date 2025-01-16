# snakefile to download data
configfile: "nbia.yaml"

COLLECTION_NAMES = config["datasets"].keys()
# COLLECTION_NAMES = list(COLLECTION_NAMES)[:1]

rule all:
  input:
    expand("procdata/{collection}.tar.gz", collection=COLLECTION_NAMES)

rule compress:
  input:
    collection_dir = "rawdata/{collection}"
  output:
    gzip_dir = "procdata/{collection}.tar.gz"
  shell:
    """
    tar -cf - {input.collection_dir} | pigz -c > {output.gzip_dir}
    """

rule download_collection:
  output:
    collection_dir = directory("rawdata/{collection}")
  script:
    "workflow/scripts/download_collection.py"