# Medical Image Test Data

This repository aims to source medical images for BHKLAB
radiomics test data.

[Releases](https://github.com/bhklab/med-image_test-data/releases)

See [this issue on Med-ImageTools](https://github.com/bhklab/med-imagetools/issues/173)
for the discussion.

Ideally, the results of this pipeline should be version controlled
and set up as release artifacts on GitHub.

## Development

### Installation

```bash
pixi install
```

### Usage

This project is just a snakemake pipeline that does the following:

1. Check the "nbia.yaml" for the list of collectiond
2. Create a job for each collection that then checks the `nbia.yaml` for the list of patients
   1. For each patient, get ALL the series
   2. Download all the series into a 'rawdata' directory
3. tar.gz each collection in the 'rawdata' directory into the `procdata` directory

```bash
pixi run snake
```

### Release

Upon a new tag being generated on github, the pipeline will be run and the results will be
packaged into a release artifact.
See [Releases](https://github.com/bhklab/med-image_test-data/releases) for all the releases.
