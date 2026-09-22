# CLAUDE.md

Guidance for Claude Code working in this repository.

## What this is

Analysis code for a continuous video–EEG rodent seizure cohort. It is **not** the modelling
repo — the training code, checkpoints and paper live in the separate `EEG/` project
(`github.com/alloydas/EEG-seizure-classification`). This repository holds the *analysis* that
sits on top: verifying the epoch-feature export, characterising the frequency bands, testing
the pre-ictal window, and ensembling stored model posteriors.

**No data is in here and none should ever be committed.** The DSI epoch export is ~7 GB and
the clip corpus ~36 GB. `.gitignore` blocks `*.edf`, `*.mp4`, `*.xlsx`, `*.npz`, `raw_epoch_features/`,
`cache/` and `output/`. Check `git status` before any `git add -A`.

## Environment

The analysis box has **numpy and scipy but no `mne`, `edfio`, `pandas` or `sklearn`** in the
base interpreter, and `environment.yml` in the sibling project warns that base is shared with
an Isaac Lab stack — **do not pip install into it**. Use a venv with `--system-site-packages`
if you need anything more. `edf/edflib.py` exists precisely because `mne` is unavailable;
prefer it over adding a dependency.

Every script reads its data root from `EEG_ROOT`, defaulting to the original box:

```bash
export EEG_ROOT=/work/mech-ai-scratch/alloy/EEG
python edf/provenance.py --animal RN244 --clips-per-day 1
```

## Layout

- `edf/` — dependency-free EDF reader, export provenance verifier, zip integrity check.
- `preictal/` — lead-seizure cohort, hour-of-day-matched controls, horizon-wise effect sizes,
  circadian permutation null, post-ictal recovery.
- `method/` — the briefs each analysis was run against. **Read `method/BRIEF.md` before
  touching the export**; it lists the gotchas below in full.
- `docs/band-structure.html` — the standalone results report.

## Data contract for the epoch export

These rules have each silently corrupted an analysis in this project. Apply all of them:

- Filter `channel_index == 1` **and** `features_valid == 1`; dedupe on
  `(recording_id, channel_index, epoch_index)` — the two-column key is one short, because one
  file exports three channels and the third is a 5 Hz high-passed copy.
- Exclude `raw_epoch_features/RoomD/RN245/` (the lab said discard it) and the `.interrupted`
  file (a byte-exact prefix of its complete twin).
- **Never `awk -F,`** — 22 files carry commas inside a quoted `source_file`, breaking 756,079
  rows. Use a real CSV parser.
- `annotation_labels_json` is a per-epoch deduplicated *set* and is **not** index-parallel to
  `annotation_ids_json`. Resolve each event's stage only from epochs where it is the sole id.
- Label an epoch ictal at `annotation_overlap_fraction >= 0.5`, never `> 0`.
- `source_unit` says `V`; the feature columns are microvolts. `supervised_label` is a constant
  placeholder, not a label.
- `features_valid = 0` blanks all 32 feature columns, and the blanking is **diurnal and
  severity-correlated** — it is a selection channel, not noise. All Spikes epochs are blanked.

## Statistical conventions — non-negotiable in this project

- **Inference is clustered on ANIMAL, never on epoch or clip.** There are ~12 M epochs; an
  epoch-level p-value is meaningless.
- **Severity claims must be reported within session as well as pooled.** Session identity
  alone — containing zero physiology — scores 0.927 pooled and exactly 0.500 within session.
  A pooled severity number is measuring which session a clip came from.
- Correlations are computed **within animal** then summarised across animals. Between-animal
  amplitude spans ~4× with no calibration, so pooling first manufactures agreement.
- The five relative powers sum to 1 and are **compositional** — use log-ratio geometry, not
  ordinary correlation.
- Band powers are heavy-tailed; work on log10.
- Prefer **level and slope** (two numbers) over five band features: the bands carry ~1.83
  effective dimensions, so five correlated tests is a multiplicity problem you do not need.
- Report counts beside every rate. Stage 5 is 14 clips on the modality intersection and 34 in
  the video validation set — do not over-read it.
- A clean null is a result. Several plausible effects here have been refuted by measurement;
  say so plainly rather than hunting for a positive.

## Established results — do not re-derive, challenge only with evidence

- The export is genuine: recomputed features match to ~1e-5, which is the clip files' 16-bit
  requantisation floor, not export error.
- Five bands ≈ two numbers (PC1 71.4%, participation ratio 1.83); a 1/f level+slope generator
  with no band coupling reproduces the whole correlation structure.
- Bands detect seizures (gamma +2.90 interictal SD) but do not grade them (within-session
  severe-vs-mild AUROC 0.46–0.55 for everything tried).
- No detectable pre-ictal change at any horizon from −60 to −1 min; powered to ~0.35 SD.
- Facility light cycle is 12:12, lights on 06:00 / off ~17:50, established three independent
  ways. `epoch_start_native` is a fixed-offset clock — do **not** localise it with a DST-aware
  timezone.

Model ensembling is a separate repository (`video-eeg-ensembling`); do not re-add it here.
