# Brief: verifying the `raw_epoch_features/` powerband export

## The claim under test
The lab sent `RoomC.zip` + `RoomD.zip` (extracted to `/work/mech-ai-scratch/alloy/EEG/raw_epoch_features/`)
with the note: *"these powerbands were pulled from the same .edf files that I provided Shreyan.
As for what each column represents are in the file attached."*

We must check (a) whether the numbers really come from those EDFs, and (b) everything a user would
need to know before trusting the file. **We never received the column-description attachment**, so the
data dictionary has to be derived empirically.

## What is already ESTABLISHED (do not re-derive; do challenge if you find contrary evidence)
* Export layout: `raw_epoch_features/{RoomC,RoomD}/[DONE_]<ANIMAL>/<day>.csv`, 1,083 CSVs, 7.1 GB,
  69 columns, one row per 10-s epoch of ONE biopotential channel.
* No full-day source EDF exists anywhere on this filesystem. The only local copy of that signal is the
  cut clips: `data_full/Data_<ANIMAL>_cropped/<day>/<clip>/eeg.edf` (+ `info.txt`).
* **Alignment key:** a clip EDF's header start time is the *source recording's* start (mne's `raw.crop`
  leaves `meas_date` alone), and `info.txt` gives the clip's wall-clock span. So export epoch at
  `epoch_start_s` maps to clip offset `epoch_start_s - (clip_start - source_start)`.
* Clip EDFs were requantized to 16 bit on write, so an exact match floors at ~1e-5 relative, never 0.
* **Feature recipe (reverse-engineered, reproduces all six power columns to 2.7e-5):**
  Welch PSD, `nperseg=4000` (4 s), Hann, 50% overlap, `detrend='constant'`, `scaling='density'`;
  band power = `np.trapezoid` over `[lo, hi]` inclusive. Bands: delta 0.5-4, theta 4-8, alpha 8-13,
  beta 13-30, gamma 30-80 Hz. `line_length_uV_per_s` = `sum(|diff|) / duration_s`.
* The export's channel is **0-based index 1 = the RAW `ECG`-labelled biopotential**, NOT the
  `ECG [FIR-HP: 5Hz]` filtered copy sitting next to it. The lab has confirmed the cohort carries no
  real ECG - it is EEG - and `fix_channel_label.py` already rewrote `channel_original_label` ECG->EEG
  on 10,992,697 rows / 979 files. `channel_identity_basis` preserves what the header actually said.
* `source_unit` is `V` in the export while the feature columns are named `_uV`; our clip EDFs are
  already in uV and match the CSV directly, so a x1e6 conversion was applied.
* Verified on RN244 (23 days, 1 clip/day): 21 days match, worst spectral rel err ~1e-3, native-clock
  delta exactly 0.0 s, `annotation_overlap_fraction` reproduces our own `info.txt` seizure boundaries
  to <5e-7. Two days yielded no fully-contained epoch; one Spikes day had every feature column blank
  (`features_valid=0`).

## Known gotchas that will silently corrupt an analysis
* 32 files have commas inside a quoted `source_file` (71 naive fields, not 69). **Never `awk -F,`** -
  use Python's `csv` module.
* `annotation_labels_json` is **not** index-parallel to `annotation_ids_json`.
* Label an epoch on `annotation_overlap_fraction >= 0.5`, never `> 0`.
* `(2)` files are **continuations, not duplicates** - dropping them costs ~127,000 epochs.
* Dedupe on `(recording_id, epoch_index)` or RN203 inflates by 47.3 h.
* `features_valid=0` rows have **blank** feature columns and are severity-correlated
  (2.97% at S2 -> 10.35% at S5 -> 100% of Spikes), so filtering on it biases the severe class.
* Exclude the one `.interrupted-*` file and the 10 header-only files from totals.
* Relative error on `mean_uV` / `skewness` / `excess_kurtosis` is meaningless when the reference value
  is near zero - normalise by `std_uV` (or by the column's own scale) before calling something a
  mismatch.

## Tools you have (already written and smoke-tested)
`<scratch>
* `edflib.py`   - dependency-free EDF/EDF+ reader (`read_header`, `read_signal`). This box has
  numpy 2.3.4 + scipy 1.16.3 but **no mne, no edfio, no pandas**. Use `python` (3.13).
* `provenance.py` - `check_clip(clip_dir, animal, day)` does the full per-clip comparison and returns a
  dict with `max_rel_exact`, `max_rel_spectral`, `max_native_clock_delta_s`,
  `max_annotation_overlap_err`, `source_file`, `source_patient`, `sample_rate_hz`, `source_unit`.
  CLI: `python provenance.py --animal RN244 [--day D] [--clips-per-day 1] [--max-clips N]` -> JSON.
  Import it (`sys.path.insert(0, ...)`) rather than rewriting it.

Write any scratch files under that same `edfcheck/` directory. Do not modify the repo.

## How to report
Return findings as structured data. A "finding" is something that would change how someone uses this
export. Distinguish hard failures (the data is wrong) from caveats (the data is fine but surprising).
Quote the exact numbers and file paths you measured - never estimate or extrapolate silently, and say
explicitly what you sampled versus what you checked exhaustively.
