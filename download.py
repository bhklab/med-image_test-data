# %%
import os
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv
from nbiatoolkit.nbia import NBIAClient
from rich import print
from yaml import safe_load

load_dotenv()

NBIA_USERNAME = os.getenv("NBIA_USERNAME")
NBIA_PASSWORD = os.getenv("NBIA_PASSWORD")


# %%
yaml_file = Path().cwd() / "nbia.yaml"
RAWDATA_PATH = Path().cwd() / "rawdata"
with yaml_file.open() as f:
    config = safe_load(f)


# %%
def download_multiple_series():
    client = NBIAClient(NBIA_USERNAME, NBIA_PASSWORD)
    for collection in config["datasets"]:
        all_series = []
        for patient in config["datasets"][collection]["patients"]:
            params = {"PatientID": patient}
            result = client.getSeries(params)
            all_series.extend(result)

        for series in all_series:
            series_uid = series["SeriesInstanceUID"]
            collection_id = series["Collection"]
            patient_id = series["PatientID"]

            zip_file = RAWDATA_PATH / collection_id / patient_id / f"{series_uid}.zip"

            csv_file = zip_file.with_suffix(".csv")
            zip_file.parent.mkdir(parents=True, exist_ok=True)
            if zip_file.exists() and csv_file.exists():
                print(f"Skipping {zip_file.name}")
                continue

            zip, metadata = client.downloadSeriesWithMetadata(
                params={"SeriesInstanceUID": series_uid},
            )

            with open(zip_file, "wb") as f:
                f.write(zip)

            # metadata is a list of dictionaries
            # convert to pandas dataframe and save
            metadata_df = pd.DataFrame(metadata)
            metadata_df.to_csv(csv_file, index=False)


download_multiple_series()
# %%
