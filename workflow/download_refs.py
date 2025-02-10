# %%
import asyncio
import re
from collections import defaultdict
from pathlib import Path
from ast import literal_eval
import pandas as pd
from nbiatoolkit import NBIA_ENDPOINT
from nbiatoolkit.nbia import NBIAClient
from nbiatoolkit.settings import Settings
from rich import print
settings = Settings()
client = NBIAClient(
    username=settings.login.nbia_username,
    password=settings.login.nbia_password,
)
import json

import numpy as np

def load_json():
  output_path = "collection_download_series.json"
  if Path(output_path).exists():

    with open(output_path, "r") as f:
        col_dict = json.load(f)
  else:
    rootdir = Path.cwd()
    csvpath = rootdir / "rtstructs_with_ref_modality.csv"
    assert csvpath.exists(), f'File not found: {csvpath}'
    df = pd.read_csv(csvpath)
    df.roi_names = df.roi_names.apply(literal_eval)

    df.head()

    to_download = df[df["ReferenceSeriesModality"].isin(["CT", "MR"])]
    to_download.head()

    to_download.collection.value_counts()

    np.random.seed(42)  # Set the random seed for reproducibility

    groups = to_download.groupby("collection")
    num_cases_to_download = 10

    col_dict = {}

    for collection, group in groups:
        
        modalities = group.ReferenceSeriesModality.unique()

        series_of_interest = []

        if len(modalities) > 1:
            MR_group = group[group.ReferenceSeriesModality == "MR"]
            CT_group = group[group.ReferenceSeriesModality == "CT"]

            # get num_cases_to_download from each modality
            series_of_interest.extend(MR_group.sample(num_cases_to_download, random_state=42).seriesuid)
            series_of_interest.extend(MR_group.sample(num_cases_to_download, random_state=42).ref_series)
            series_of_interest.extend(CT_group.sample(num_cases_to_download, random_state=42).seriesuid)
            series_of_interest.extend(CT_group.sample(num_cases_to_download, random_state=42).ref_series)

        else:
            series_of_interest.extend(group.sample(num_cases_to_download, random_state=42).seriesuid)
            series_of_interest.extend(group.sample(num_cases_to_download, random_state=42).ref_series)
        
        col_dict[collection] = series_of_interest

    with open(output_path, "w") as f:
        json.dump(col_dict, f, indent=4)



  return col_dict


async def download_and_save_series(
  client: NBIAClient, SeriesInstanceUID: str, zip_file: Path
):
  data_task = client._downloadSeries(params={"SeriesInstanceUID": SeriesInstanceUID})
  bytes_data = await data_task
  with open(zip_file, "wb") as f:
    f.write(bytes_data)
  return SeriesInstanceUID


async def main():
  col_dict = load_json()

  for collection, series_list in col_dict.items():
    OUTPUT_DIR = Path(f"med_image_net/{collection}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    tasks = []
    for series in series_list:
      zip_file = OUTPUT_DIR / f"{series}.zip"
      if zip_file.exists():
        continue
      tasks.append(download_and_save_series(client, series, zip_file))

    await asyncio.gather(*tasks)

if __name__ == "__main__":
  asyncio.run(main())