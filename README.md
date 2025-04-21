# Medical Image Test Data

This repository aims to source medical images for BHKLAB
radiomics test data.

[Public Releases](https://github.com/bhklab/med-image_test-data/releases)
[Private Releases](https://github.com/bhklab/med-image_test-data_private/releases)

See [this issue on Med-ImageTools](https://github.com/bhklab/med-imagetools/issues/173)
for the discussion.

## Development

### Installation

```bash
pixi install
```

### Usage

This project is just a snakemake pipeline that does the following:

1. Check the "nbia_datasets.yaml" for the list of collectiond
2. Create a job for each collection that then checks the `nbia.yaml` for the list of patients
   1. For each patient, get ALL the series
   2. Download all the series into a 'rawdata' directory
3. tar.gz each collection in the 'rawdata' directory into the `procdata` directory

### Access Type 

This pipeline is by default, using the 'public' access type.
The bhklab has a clone of this repo that runs the pipeline in `private` mode.
To run the pipeline in private mode, you need to set up your environment variables:

export the username and password AND create a `.env` file with the credentials.

```bash
export NBIA_USERNAME=""
export NBIA_PASSWORD=""
```

Create a `.env` file like:
```
LOGIN__NBIA_USERNAME=""
LOGIN__NBIA_PASSWORD=""
```

You also need to set the `ACCESS_TYPE` to `private` 

```console
export ACCESS_TYPE=private
```

### Running the pipeline
Then:

```bash
pixi run snake
```

### Release

Upon a new tag being generated on github, the pipeline will be run and the results will be
packaged into a release artifact.
See [Releases](https://github.com/bhklab/med-image_test-data/releases) for all the releases.
