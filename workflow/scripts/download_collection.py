import os
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv
from nbiatoolkit.nbia import NBIAClient
from rich import print
import asyncio

load_dotenv()


NBIA_USERNAME = os.getenv("NBIA_USERNAME")
NBIA_PASSWORD = os.getenv("NBIA_PASSWORD")
client = NBIAClient(NBIA_USERNAME, NBIA_PASSWORD)


COLLECTION = snakemake.wildcards.collection

config = snakemake.config["datasets"][COLLECTION]

print(config)

print(f"Downloading data for {COLLECTION}")

OUTPUT_DIR = Path(snakemake.output[0])
series_list = []

for patient in config["patients"]:
    params = {"PatientID": patient}
    result = client.getSeries(params)
    series_list.extend(result)
print(f"Found {len(series_list)} series")

async def download_multiple_series(SeriesInstanceUID: str) -> dict[str, bytes]:
  data_task = client._downloadSeries(
    params={"SeriesInstanceUID": SeriesInstanceUID}
  )
  bytes_data = await data_task
  return {SeriesInstanceUID: bytes_data}

async def main():
  tasks = []
  for series in series_list:
    zip_file = OUTPUT_DIR / series["PatientID"] / f"{series['SeriesInstanceUID']}.zip"
    if zip_file.exists():
      print(f"Skipping {zip_file.name}")
      continue
    tasks.append(download_multiple_series(series["SeriesInstanceUID"]))

  results = await asyncio.gather(*tasks)

  for result in results:
    for SeriesInstanceUID, zip_data in result.items():
      series = next(s for s in series_list if s["SeriesInstanceUID"] == SeriesInstanceUID)
      zip_file = OUTPUT_DIR / series["PatientID"] / f"{SeriesInstanceUID}.zip"
      zip_file.parent.mkdir(parents=True, exist_ok=True)
      with open(zip_file, "wb") as f:
        f.write(zip_data)

asyncio.run(main())
