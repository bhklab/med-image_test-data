# %%
import asyncio
from io import BytesIO
import re
from collections import defaultdict
from itertools import chain
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import pandas as pd
from nbiatoolkit import NBIA_ENDPOINT
from nbiatoolkit import nbia
from nbiatoolkit.settings import Settings
from nbiatoolkit import dicomtags as dcm
import pydicom
from rich import print


# %%
settings = Settings()
client = nbia.NBIAClient(
    username=settings.login.nbia_username,
    password=settings.login.nbia_password,
)
try:
    loop = asyncio.get_running_loop()
    result = loop.run_until_complete(client.query(NBIA_ENDPOINT.GET_COLLECTIONS))
except RuntimeError:  # no event loop running
    result = asyncio.run(client.query(NBIA_ENDPOINT.GET_COLLECTIONS))

# %%
collection_list = [
    "C4KC-KiTS",
    "CPTAC-CCRCC",
    "CPTAC-HNSCC",
    "CPTAC-PDA",
    "CPTAC-SAR",
    "CPTAC-UCEC",
    "CT Lymph Nodes",
    "HEAD-NECK-RADIOMICS-HN1",
    "HNSCC",
    "HNSCC-3DCT-RT",
    "Head-Neck-PET-CT",
    "Pancreas-CT",
    "Pancreatic-CT-CBCT-SEG",
    "QIN-HEADNECK",
    "RADCURE",
    "TCGA-BLCA",
    "TCGA-HNSC",
    "TCGA-KICH",
    "TCGA-KIRC",
    "TCGA-KIRP",
    "TCGA-LIHC",
    "TCGA-OV",
    "TCGA-STAD",
]

all_results = defaultdict(list)

all_series_path = Path("all_series.csv")

if all_series_path.exists():
    all_series_df = pd.read_csv(all_series_path)
else:
    async def fetch_series():
        tasks = [
            client._getSeries(params={"Collection": collection})
            for collection in collection_list
        ]
        
        results = await asyncio.gather(*tasks)
        return list(chain(*results))

    all_series = asyncio.run(fetch_series())

    all_series_df = pd.DataFrame(all_series)
    columns_of_interest = [
        "Collection",
        "PatientID",
        "StudyInstanceUID",
        "Modality",
        "SeriesInstanceUID",
        "SeriesNumber",
        "BodyPartExamined",
        "ImageCount",
        "TimeStamp",
        "CollectionURI",
        "FileSize",
        "DateReleased",
    ]

    all_series_df = all_series_df[columns_of_interest]
    all_series_df.reset_index(drop=True, inplace=True)
    all_series_df.to_csv(all_series_path, index=False)

modalities_of_interest = ["CT", "MR", "RTSTRUCT", "SEG"]
series_df = all_series_df[all_series_df["Modality"].isin(modalities_of_interest)]

print(pd.crosstab(all_series_df["Collection"], all_series_df["Modality"]).sort_values(["RTSTRUCT"], ascending=False, inplace=False))
# %%

print(pd.crosstab(series_df["Collection"], series_df["Modality"]).sort_values(["RTSTRUCT"], ascending=False, inplace=False))
# %%
groups = series_df.groupby(by="Collection")


col = "Head-Neck-PET-CT"

# get group:
group = groups.get_group(col)

series_list = group[group["Modality"] == "RTSTRUCT"]['SeriesInstanceUID'].to_list()
print(f"Number of RTSTRUCT series in {col}: {len(series_list)}")


def load_dicom_series_from_zip(zip_data: bytes) -> pydicom.dataset.FileDataset:
    """
    Load a DICOM series from a ZIP archive in memory and return a SimpleITK image.

    Parameters
    ----------
    zip_data : bytes
        The content of the ZIP archive containing the DICOM series.

    Returns
    -------
    sitk.Image
        A SimpleITK image created from the DICOM series in the ZIP archive.
    """
    with TemporaryDirectory() as temporary_dir:
        temp_dir = Path(temporary_dir)

        try:
            with ZipFile(BytesIO(zip_data)) as zf:
                dicom_files = [name for name in zf.namelist() if name.endswith(".dcm")]
                for dicom_file in dicom_files:
                    with zf.open(dicom_file) as file:
                        temp_file_path = temp_dir / Path(dicom_file).name
                        with temp_file_path.open("wb") as temp_file:
                            temp_file.write(file.read())
            
            dicom_files = list(temp_dir.glob("*.dcm"))
            assert dicom_files, "No DICOM files found in the ZIP archive."
            assert len(dicom_files) == 1, f"Only one DICOM file per series is supported. Found {len(dicom_files)} DICOM files."
            dicom_file = dicom_files[0]
            dicom = pydicom.dcmread(dicom_file)
        finally:
            shutil.rmtree(temp_dir)
        return dicom

    return image

# data = asyncio.run(client._downloadSeries(params = {"SeriesInstanceUID": series_list[0]}))
def download_RTSTRUCTS():
    for group_name, group in groups:

        async def download_series(row: pd.Series):
            file_path = Path("rawdata/rtstructs/") / row.Collection / f"{row.SeriesInstanceUID}.dcm"
            if file_path.exists():
                return file_path
            file_path.parent.mkdir(parents=True, exist_ok=True)
            data = await client._downloadSeries(params={"SeriesInstanceUID": row.SeriesInstanceUID})
            dicom = load_dicom_series_from_zip(data)
            dicom.save_as(file_path)

            return file_path


        rtstructs = group[group["Modality"] == "RTSTRUCT"]

        if rtstructs.empty:
            client.progress_bar.console.print(f"[bold red]No RTSTRUCT series found in {group_name}")
            continue

        async def download_batch(batch):
            tasks = [download_series(row) for idx, row in batch.iterrows()]
            return await asyncio.gather(*tasks)

        batch_size = 25
        task = client.progress_bar.add_task(f"Downloading RTSTRUCTs from {group_name}", total=len(rtstructs))
        for i in range(0, len(rtstructs), batch_size):
            batch = rtstructs.iloc[i:i + batch_size]
            paths = asyncio.run(download_batch(batch))
            client.progress_bar.update(task, advance=len(paths))


def download_SEG():
    for group_name, group in groups:
        async def download_series(row: pd.Series):
            file_path = Path("rawdata/seg/") / row.Collection / f"{row.SeriesInstanceUID}.dcm"
            if file_path.exists():
                return file_path
            file_path.parent.mkdir(parents=True, exist_ok=True)
            data = await client._downloadSeries(params={"SeriesInstanceUID": row.SeriesInstanceUID})
            dicom = load_dicom_series_from_zip(data)
            dicom.save_as(file_path)

            return file_path


        seg_series = group[group["Modality"] == "SEG"]

        if seg_series.empty:
            client.progress_bar.console.print(f"[bold red]No SEG series found in {group_name}")
            continue

        async def download_batch(batch):
            tasks = [download_series(row) for idx, row in batch.iterrows()]
            return await asyncio.gather(*tasks)
        
        batch_size = 25
        task = client.progress_bar.add_task(f"Downloading SEGs from {group_name}", total=len(seg_series))
        for i in range(0, len(seg_series), batch_size):
            batch = seg_series.iloc[i:i + batch_size]
            paths = asyncio.run(download_batch(batch))
            client.progress_bar.update(task, advance=len(paths))

if __name__ == "__main__":
    # download_RTSTRUCTS()
    download_SEG()

# dicom = load_dicom_series_from_zip(data)



# # %%
# result= nbia.get_single_image_per_series(
#     collection="NSCLC-Radiomics",
#     modality="RTSTRUCT",
#     num_series=1,
#     client=client,
# )
# client.progress_bar.console.print(result)