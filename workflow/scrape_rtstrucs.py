from imgtools.dicom import load_rtstruct_dcm, extract_roi_names, rtstruct_reference_uids
from imgtools.dicom.input.dicom_reader import DicomInput
from imgtools import find_dicoms
from dataclasses import dataclass, asdict
from pathlib import Path
@dataclass(frozen=True, slots=True)
class RTStruct:
  roi_names: list[str]
  ref_series: str
  ref_study: str

  @classmethod
  def from_rtstruct(cls, rtstruct: DicomInput):
    dcm = load_rtstruct_dcm(rtstruct)
    roi_names = extract_roi_names(dcm)
    ref_series, ref_study = rtstruct_reference_uids(dcm)
    return cls(roi_names, ref_series, ref_study)
  
  def to_dict(self):
    return asdict(self)


def get_rtstructs(rtstructs: list[DicomInput]) -> dict[str, RTStruct]:
  # return [RTStruct.from_rtstruct(rtstruct) for rtstruct in rtstructs]
  rtstructs = {rtstruct: RTStruct.from_rtstruct(rtstruct) for rtstruct in rtstructs}
  return rtstructs


def get_rtstructs_from_dir(rtstruct_dir: str) -> dict[str, RTStruct]:
  rtdur = Path(rtstruct_dir)

  assert rtdur.exists() and rtdur.is_dir(), f"{rtdur} does not exist or is not a directory"
  rtstructs = find_dicoms(rtdur, recursive=True, check_header=False)

  return get_rtstructs(rtstructs)

if __name__ == "__main__":
  import pandas as pd 
  from tqdm import tqdm
  from rich.progress import track
  COLLECTION_LIST = Path("rawdata/rtstructs").glob("*")

  OUTPUT_PATH = "rtstructs.csv"

  Path(OUTPUT_PATH).unlink(missing_ok=True)

  for COLLECTION in track(COLLECTION_LIST, description="Processing collections..."):
    rtstructs = {
      Path(series_path).stem: rt.to_dict()
      for series_path, rt in get_rtstructs_from_dir(f"rawdata/rtstructs/{COLLECTION.name}").items()
    }

    # create a row for each seriesuid and the corresponding roi_names, ref_series, ref_study
    df = pd.DataFrame(rtstructs).T.reset_index()
    df.columns = ["seriesuid", "roi_names", "ref_series", "ref_study"]
    df["collection"] = COLLECTION.name
    # make collection the first col
    cols = ["collection", "seriesuid", "ref_series", "ref_study", "roi_names"]
    df = df[cols]

    df.to_csv(OUTPUT_PATH, mode='a', header=not Path(OUTPUT_PATH).exists(), index=False)


