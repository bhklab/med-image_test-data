# snakefile to download data
configfile: "nbia.yaml"

COLLECTION_NAMES = config["datasets"].keys()
# COLLECTION_NAMES = list(COLLECTION_NAMES)[:1]

rule all:
  input:
    expand("results/{collection}.tar.gz", collection=COLLECTION_NAMES),
    expand("results/{collection}_summary.md", collection=COLLECTION_NAMES)
  shell:
    "echo 'All done!'"

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
    summary_file = "results/{collection}_summary.md"
  script:
    "workflow/scripts/summarize_metadata.py"

rule unzip:
  input:
    collection_dir = "rawdata/{collection}"
  output:
    unzip_dir = directory("procdata/{collection}")
  log:
    "logs/unzip_{collection}.log"
  shell:
    # every file in the rawdata directory is a .zip of the series 
    # rawdata / {collection} / {a patient id} / {many series id}.zip
    # iterate through all the zip files and unzip them
    # while maintaining the directory structure
    # into the procdata directory
    # exec > >(tee -a logfile.log) 2> >(tee -a logfile.log >&2)
    """
    for zip_file in $(find {input.collection_dir} -name '*.zip'); do
      RELATIVE_DIR=$(basename $(dirname "$zip_file"))
      FILENAME=$(basename "$zip_file" .zip)
      TARGET_DIR={output.unzip_dir}/"$RELATIVE_DIR"/"$FILENAME"
      mkdir -p "$TARGET_DIR"
      unzip -o "$zip_file" -d "$TARGET_DIR" >> {log} 2>&1
      find "$TARGET_DIR" -name 'LICENSE' -delete
    done
    """

rule download_collection:
  output:
    collection_dir = directory("rawdata/{collection}"),
    metadata_file = "metadata/{collection}.json"
  script:
    "workflow/scripts/download_collection.py"