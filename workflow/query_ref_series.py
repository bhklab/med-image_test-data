import pandas as pd
from ast import literal_eval
from tqdm import tqdm

all_series = pd.read_csv("all_series.csv")

cols = ["Collection", "Modality", "SeriesInstanceUID"]

all_series = all_series[cols]

print(all_series["Collection"].value_counts())


rtstructs = pd.read_csv("rtstructs.csv")
rtstructs.roi_names = rtstructs.roi_names.apply(literal_eval)

series_list = rtstructs["seriesuid"].to_list()

unique_s = set(series_list)

print(f"Number of unique series: {len(unique_s)} out of {len(series_list)}")


references = rtstructs.ref_series.to_list()

unique_r = set(references)

print(f"Number of unique references: {len(unique_r)} out of {len(references)}")


all_series_uids = all_series["SeriesInstanceUID"].to_list()

unique_all = set(all_series_uids)

print(f"Number of unique all series: {len(unique_all)} out of {len(all_series_uids)}")

print(f"Number of unique series in references: {len(unique_r.intersection(unique_all))}")
print(f"Missing series in references: {len(unique_r.difference(unique_all))}")

# get the row in all_series that matches each reference, then add a column to rtstructs with the 'row.Modality' value

ref_mods = []

for ref in tqdm(unique_r.intersection(unique_all)):
    row = all_series[all_series["SeriesInstanceUID"] == ref]
    row_dict = row.to_dict("records")[0]

    ref_mod = row_dict["Modality"]
    rtstructs.loc[rtstructs["ref_series"] == ref, "ReferenceSeriesModality"] = ref_mod
    ref_mods.append(ref_mod)
    # break


# get counts of each modality in the references
ref_mods = pd.Series(ref_mods)
print(ref_mods.value_counts())


rtstructs.to_csv("rtstructs_with_ref_modality.csv", index=False)
# to all_series, add 