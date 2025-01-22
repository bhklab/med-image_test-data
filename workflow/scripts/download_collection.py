import asyncio
import json
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

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
            print(f"Skipping {zip_file.name}")
            continue
        tasks.append(download_multiple_series(client, series["SeriesInstanceUID"]))

    results = await asyncio.gather(*tasks)

    for result in results:
        for SeriesInstanceUID, zip_data in result.items():
            series = next(
                s for s in series_list if s["SeriesInstanceUID"] == SeriesInstanceUID
            )
            zip_file = OUTPUT_DIR / series["PatientID"] / f"{SeriesInstanceUID}.zip"
            zip_file.parent.mkdir(parents=True, exist_ok=True)
            with open(zip_file, "wb") as f:
                f.write(zip_data)


if __name__ == "__main__":
    COLLECTION = snakemake.wildcards.collection
    config = snakemake.config["datasets"][COLLECTION]
    METADATA_FILE = Path(snakemake.output["metadata_file"])
    OUTPUT_DIR = Path(snakemake.output["collection_dir"])

    # Your script logic here
    settings = Settings()
    client = NBIAClient(
        username=settings.login.nbia_username,
        password=settings.login.nbia_password,
    )

    print(config)

    print(f"Downloading data for {COLLECTION}")

    series_list = []

    for patient in config["patients"]:
        params = {"PatientID": patient}
        result = client.getSeries(params)
        series_list.extend(result)
    print(f"Found {len(series_list)} series")

    json_file_content = json.dumps(series_list, indent=4)
    with METADATA_FILE.open("w") as f:
        f.write(json_file_content)

    asyncio.run(main(client, series_list, OUTPUT_DIR))
