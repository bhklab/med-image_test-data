# %%
import asyncio
import re
from collections import defaultdict
from itertools import chain
from pathlib import Path

import pandas as pd
from nbiatoolkit import NBIA_ENDPOINT
from nbiatoolkit.nbia import NBIAClient
from nbiatoolkit.settings import Settings
from nbiatoolkit import dicomtags as dcm
from rich import print

# %%
settings = Settings()
client = NBIAClient(
    username=settings.login.nbia_username,
    password=settings.login.nbia_password,
)

client.query(NBIA_ENDPOINT.GET_COLLECTIONS)
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

# %%

pd.crosstab(all_series_df["Collection"], all_series_df["Modality"]).sort_values(["RTSTRUCT"], ascending=False, inplace=False)
# %%

pd.crosstab(series_df["Collection"], series_df["Modality"]).sort_values(["RTSTRUCT"], ascending=False, inplace=False)

# %%


rtstructs = series_df[series_df["Modality"] == "RTSTRUCT"]
rtstructs["Collection"].value_counts()
rtstructs_backup = rtstructs.copy()



# %%

# # Function to get ReferencedSeriesUID
# async def get_referenced_series_uid(client, series_uid):
#     # Query the DICOM tags
#     tags = await client.query(NBIA_ENDPOINT.GET_DICOM_TAGS, {"SeriesUID": series_uid})
#     tags_df = pd.DataFrame(tags)

#     # Extract the ReferencedSeriesUID
#     return dcm.tags.getReferencedSeriesUIDS(tags_df)

# # Wrapper function to check and apply the ReferencedSeriesUID
# async def add_referenced_series_uid(row, client):
#     if pd.notnull(row.get("ReferencedSeriesUID")):  # Skip if already exists
#         return row["ReferencedSeriesUID"]
    
#     # Call async function to get ReferencedSeriesUID
#     return await get_referenced_series_uid(client, row["SeriesInstanceUID"])

# # Update DataFrame with new column in place
# async def add_referenced_column_in_place(client, rtstructs, max_concurrent_tasks=10, n=None):
#     # Helper function to run tasks with a limit on concurrency
#     async def limited_as_completed(tasks, limit):
#         semaphore = asyncio.Semaphore(limit)

#         async def sem_task(coro):
#             async with semaphore:
#                 return await coro

#         # Use `gather` to run tasks with concurrency control
#         return await asyncio.gather(*(sem_task(task) for task in tasks))
#     if n is not None:
#         rtstructs = rtstructs.head(n)
#     # Apply async function row-wise with concurrency limit
#     tasks = [
#         add_referenced_series_uid(row, client) for _, row in rtstructs.iterrows()
#     ]

#     subset_df = rtstructs.head(len(tasks))
#     referenced_series_uids = await limited_as_completed(tasks, max_concurrent_tasks)
#     # Update the subset of DataFrame
#     subset_df.loc[:, "ReferencedSeriesUID"] = referenced_series_uids

#     return subset_df

# # Example usage
# async def main():
#     rtstructs_backup["ReferencedSeriesUID"] = None
#     # Iterate over each collection and process separately
#     for collection in rtstructs_backup["Collection"].unique():
#         client.progress_bar.console.print(f"Processing collection: {collection}")
#         collection_rtstructs = rtstructs_backup[rtstructs_backup["Collection"] == collection]
#         refd_rtstructs = await add_referenced_column_in_place(client, collection_rtstructs, max_concurrent_tasks=10, n=10)
        
#         # Update the original DataFrame
#         rtstructs_backup.update(refd_rtstructs)


# Run the main function
asyncio.run(main())

# %%
# print(all_series_df.query("SeriesInstanceUID in @refd_rtstructs['ReferencedSeriesUID']"))
print(rtstructs_backup[rtstructs_backup.notnull().all(axis=1)].ReferencedSeriesUID)

# rt_tags_dict = defaultdict(list)
# for collection, df in rtstructs_grouped.groupby("Collection"):
#     rt_tags_dict[collection] = [result for result in results if result["SeriesInstanceUID"].iloc[0] in df["SeriesInstanceUID"].values]



# %%
