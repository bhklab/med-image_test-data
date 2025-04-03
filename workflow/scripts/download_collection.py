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


async def download_and_extract_series(
    client: NBIAClient, series: str, output_path: Path
) -> bool:
    """Download and extract a series in one operation to save memory."""
    logger.info(f"Downloading {series}")
    try:

        zip_data = await client._download_series(series)
        logger.info(f"Extracting {series}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_data) as z:
            z.extractall(output_path)
        return True
    except Exception as e:
        logger.error(f"Error processing series {series}: {e}")
        return False


async def main(client: NBIAClient, series_list: pd.DataFrame, OUTPUT_DIR: Path):
    # filter out series where the `FileSize` column is greater than 1GB
    too_big = series_list[series_list["FileSize"] > 1_000_000_000]
    if len(too_big) > 0:
        logger.warning(
            f"Found {len(too_big)} series with size greater than 1GB. Skipping them."
        )
        logger.warning(f"Too Big Series: {too_big.SeriesInstanceUID.tolist()}")
    series_list = series_list[series_list["FileSize"] <= 1_000_000_000]

    series_list = series_list.SeriesInstanceUID.unique()

    logger.info(f"Found {len(series_list)} series")

    # Create a directory for each patient
    for patient_id in metadata_df.PatientID.unique():
        patient_dir = OUTPUT_DIR / patient_id
        patient_dir.mkdir(parents=True, exist_ok=True)

    # add a column for the path to the zip file
    metadata_df["zip_file"] = metadata_df.apply(
        lambda row: OUTPUT_DIR
        / row.PatientID
        / f"{row.Modality}_Series{row.SeriesInstanceUID[-8:]}",
        axis=1,
    )
    metadata_df.set_index("SeriesInstanceUID", inplace=True)

    # Download and extract the series concurrently
    tasks = []

    for series in series_list:
        if metadata_df.loc[series].zip_file.exists():
            logger.warning(f"Skipping {series}")
            continue

        # Download and extract the series in one go
        tasks.append(
            download_and_extract_series(
                client, series, metadata_df.loc[series].zip_file
            )
        )


    results = await asyncio.gather(*tasks)

    for result in results:
        if result is False:
            logger.error("Failed to download series something ")
            import sys

            sys.exit(1)
    successful_downloads = sum(results)
    logger.info(
        f"Successfully processed {successful_downloads} out of {len(tasks)} series"
    )


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
    client.disable_progress_bar = True  # doesnt work right now lol

    # print(f"Downloading data for {COLLECTION}")
    logger.info(f"Downloading data for {COLLECTION}")

    # Read the metadata file
    metadata_df = pd.read_csv(METADATA_FILE)

    # json_file_content = json.dumps(series_list, indent=4)
    # with METADATA_FILE.open("w") as f:
    #     f.write(json_file_content)

    asyncio.run(main(client, metadata_df, OUTPUT_DIR))
