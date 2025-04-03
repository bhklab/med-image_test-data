import asyncio
import io
import json
import logging
import warnings
import zipfile
from pathlib import Path
from typing import TYPE_CHECKING

from nbiatoolkit.logging_config import logger
from nbiatoolkit.nbia import NBIAClient
from nbiatoolkit.settings import Settings
from nbiatoolkit.models.nbia_responses import Series, SeriesList
from rich import print

if TYPE_CHECKING:
    from snakemake.script import snakemake  # type: ignore

# Suppress SyntaxWarning
warnings.filterwarnings("ignore", category=SyntaxWarning)

if __name__ == "__main__":
    COLLECTION = snakemake.wildcards.collection
    config = snakemake.config["datasets"][COLLECTION]
    METADATA_FILE = Path(snakemake.output["metadata_file"])
    # OUTPUT_DIR = Path(snakemake.output["collection_dir"])
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

    param_list = []

    for key in config:
        # could be PatientID, StudyInstanceUID, SeriesInstanceUID
        assert key in [
            "PatientID",
            "StudyInstanceUID",
            "SeriesInstanceUID",
        ], (
            f"Invalid key {key} in config. Expected one of PatientID, StudyInstanceUID, SeriesInstanceUID"
        )
        param_list.append({key: identifier for identifier in config[key]})

    # param_list = [
    #     {"PatientID": patient}
    #     for patient in config["patients"]
    # ]

    series = client.getSeries(param_list)

    logger.info(f"Found {len(series)} series")

    s: SeriesList = Series.from_dicts(series)

    s.df.to_csv(METADATA_FILE, index=False)
