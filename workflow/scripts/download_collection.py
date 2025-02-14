import asyncio
import io
import json
import logging
import warnings
import zipfile
from pathlib import Path
from typing import TYPE_CHECKING

from nbiatoolkit import logger
from nbiatoolkit.nbia import NBIAClient
from nbiatoolkit.settings import Settings
from rich import print

if TYPE_CHECKING:
    from snakemake.script import snakemake  # type: ignore


# Suppress SyntaxWarning
warnings.filterwarnings("ignore", category=SyntaxWarning)


async def download_multiple_series(
    client: NBIAClient, SeriesInstanceUID: str
) -> dict[str, bytes]:
    data_task = client._downloadSeries(params={"SeriesInstanceUID": SeriesInstanceUID})
    bytes_data = await data_task
    return {SeriesInstanceUID: bytes_data}


async def main(client: NBIAClient, series_list: list, OUTPUT_DIR: Path):
    tasks = []
    for series in series_list:
        zip_file = (
            OUTPUT_DIR / series["PatientID"] / f"{series['SeriesInstanceUID']}.zip"
        )
        if zip_file.exists():
            logger.warning(f"Skipping {zip_file.name}")
            continue
        tasks.append(download_multiple_series(client, series["SeriesInstanceUID"]))

    results = await asyncio.gather(*tasks)
    logger.info(f"Downloaded {len(results)} series")
    for result in results:
        for SeriesInstanceUID, zip_data in result.items():
            series = next(
                s for s in series_list if s["SeriesInstanceUID"] == SeriesInstanceUID
            )
            logger.info(f"Extracting {series['SeriesInstanceUID']}")
            series_dir = OUTPUT_DIR / series["PatientID"] / SeriesInstanceUID
            series_dir.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(io.BytesIO(zip_data)) as z:
                z.extractall(series_dir)


if __name__ == "__main__":
    COLLECTION = snakemake.wildcards.collection
    config = snakemake.config["datasets"][COLLECTION]
    METADATA_FILE = Path(snakemake.output["metadata_file"])
    OUTPUT_DIR = Path(snakemake.output["collection_dir"])
    LOG_FILE = Path(snakemake.log[0])

    # add a file handler to the logger
    file_handler = logging.FileHandler(LOG_FILE)
    logger.addHandler(file_handler)

    # Your script logic here
    settings = Settings()
    client = NBIAClient(
        username=settings.login.nbia_username,
        password=settings.login.nbia_password,
    )

    # print(config)

    # print(f"Downloading data for {COLLECTION}")
    logger.info(f"Downloading data for {COLLECTION}")

    series_list = []

    for patient in config["patients"]:
        params = {"PatientID": patient}
        result = client.getSeries(params)
        series_list.extend(result)
    logger.info(f"Found {len(series_list)} series")

    json_file_content = json.dumps(series_list, indent=4)
    with METADATA_FILE.open("w") as f:
        f.write(json_file_content)

    asyncio.run(main(client, series_list, OUTPUT_DIR))
