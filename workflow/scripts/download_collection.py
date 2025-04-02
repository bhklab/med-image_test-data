import asyncio
import io
import json
import logging
import warnings
import zipfile
from pathlib import Path
from typing import TYPE_CHECKING
import pandas as pd
from nbiatoolkit import logger
from nbiatoolkit.nbia import NBIAClient
from nbiatoolkit.settings import Settings
from rich import print

if TYPE_CHECKING:
    from snakemake.script import snakemake  # type: ignore


# Suppress SyntaxWarning
warnings.filterwarnings("ignore", category=SyntaxWarning)


async def main(client: NBIAClient, series_list: pd.DataFrame, OUTPUT_DIR: Path):
    series_list = metadata_df.SeriesInstanceUID.unique()

    logger.info(f"Found {len(series_list)} series")

    # Create a directory for each patient
    for patient_id in metadata_df.PatientID.unique():
        patient_dir = OUTPUT_DIR / patient_id
        patient_dir.mkdir(parents=True, exist_ok=True)

    # add a column for the path to the zip file
    metadata_df["zip_file"] = metadata_df.apply(
        lambda row: OUTPUT_DIR / row.PatientID / f"{row.Modality}_Series{row.SeriesInstanceUID[-8:]}.zip", axis=1
    )
    metadata_df.set_index("SeriesInstanceUID", inplace=True)

    # Download the series
    tasks = []

    async def series_map(series)-> dict[str: bytes]:
        logger.info(f"Downloading {series}")
        zip_data = await client._download_series(series)
        return {series: zip_data}

    for series in series_list:
        if metadata_df.loc[series].zip_file.exists():
            logger.warning(f"Skipping {series}")
            continue

        # download the series
        tasks.append(series_map(series))

    results = await asyncio.gather(*tasks)

    logger.info(f"Downloaded {len(results)} series")
    for result in results:
        for SeriesInstanceUID, zip_data in result.items():
            series = metadata_df.loc[SeriesInstanceUID]
            logger.info(f"Extracting {SeriesInstanceUID}")
            zip_file = series.zip_file
            zip_file.parent.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(zip_data) as z:
                z.extractall(zip_file)

if __name__ == "__main__":
    COLLECTION = snakemake.wildcards.collection
    config = snakemake.config["datasets"][COLLECTION]
    METADATA_FILE = Path(snakemake.input["metadata_file"])
    OUTPUT_DIR = Path(snakemake.output["collection_dir"])
    LOG_FILE = Path(snakemake.log[0])

    # add a file handler to the logger
    file_handler = logging.FileHandler(LOG_FILE)
    # remove existing handlers
    logger.handlers = []
    logger.addHandler(file_handler)
    # format with log level and time
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)

    # Your script logic here
    settings = Settings()

    client = NBIAClient.from_settings(settings)
    client.disable_progress_bar = True # doesnt work right now lol

    # print(f"Downloading data for {COLLECTION}")
    logger.info(f"Downloading data for {COLLECTION}")

    # Read the metadata file
    metadata_df = pd.read_csv(METADATA_FILE)

    # json_file_content = json.dumps(series_list, indent=4)
    # with METADATA_FILE.open("w") as f:
    #     f.write(json_file_content)

    asyncio.run(main(client, metadata_df, OUTPUT_DIR))
