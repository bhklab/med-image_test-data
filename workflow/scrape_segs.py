from dataclasses import dataclass, asdict
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from rich.progress import track

from imgtools.dicom import load_seg_dcm
from imgtools.dicom.input.dicom_reader import DicomInput
from imgtools import find_dicoms


def extract_segment_names(dcm):
  return [seg.SegmentLabel for seg in dcm.SegmentSequence]


def seg_reference_uids(dcm):
  return dcm.ReferencedSeriesSequence[0].SeriesInstanceUID

@dataclass(frozen=True, slots=True)
class SEGModality:
	segment_names: list[str]
	ref_series: str

	@classmethod
	def from_seg(cls, seg: DicomInput):
		dcm = load_seg_dcm(seg)
		segment_names = extract_segment_names(dcm)
		ref_series = seg_reference_uids(dcm)
		return cls(segment_names, ref_series)

	def to_dict(self):
		return asdict(self)


def get_segmodalities(segs: list[DicomInput]) -> dict[str, SEGModality]:
	segmodalities = {str(seg): SEGModality.from_seg(seg) for seg in segs}
	return segmodalities


def get_segmodalities_from_dir(seg_dir: str) -> dict[str, SEGModality]:
	seg_path = Path(seg_dir)

	assert seg_path.exists() and seg_path.is_dir(), f"{seg_path} does not exist or is not a directory"
	segs = find_dicoms(seg_path, recursive=True, check_header=False)

	return get_segmodalities(segs)


if __name__ == "__main__":
	SEG_COLLECTION_LIST = Path("rawdata/seg").glob("*")
	SEG_OUTPUT_PATH = "segmodalities.csv"

	Path(SEG_OUTPUT_PATH).unlink(missing_ok=True)

	for COLLECTION in track(SEG_COLLECTION_LIST, description="Processing collections..."):
		segmodalities = {
			Path(series_path).stem: seg.to_dict()
			for series_path, seg in get_segmodalities_from_dir(f"rawdata/seg/{COLLECTION.name}").items()
		}

		# create a row for each seriesuid and the corresponding segment names, ref_series
		df = pd.DataFrame(segmodalities).T.reset_index()
		df.columns = ["seriesuid", "segment_names", "ref_series"]
		df["collection"] = COLLECTION.name
		# make collection the first column
		cols = ["collection", "seriesuid", "ref_series","segment_names"]
		df = df[cols]

		df.to_csv(SEG_OUTPUT_PATH, mode='a', header=not Path(SEG_OUTPUT_PATH).exists(), index=False)