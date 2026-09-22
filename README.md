# Rodent EEG band analysis

Analysis code and results for a continuous video–EEG rodent seizure cohort: verifying a
33,455-hour DSI epoch-feature export, characterising how the classical frequency bands
relate to each other, and testing what they do around seizures.

Two self-contained HTML reports (open either in a browser):

- **[Five Bands, Two Numbers](docs/band-structure.html)** — how the frequency bands relate to
  each other, and what they do around seizures.
- **[Forty Runs, No Retraining](docs/ensembling.html)** — ensembling the stored model
  posteriors: free gains for video, none for EEG.

## What is here

| path | what it does |
|---|---|
| `edf/edflib.py` | Dependency-free EDF / EDF+ reader (numpy only). Written because the analysis box has numpy + scipy but no `mne` or `edfio`. |
| `edf/provenance.py` | Verifies that an epoch-feature export was computed from the EDFs it claims, by recomputing features from the source samples and comparing column by column. |
| `edf/zipcheck.py` | Byte-exact check of an extracted delivery against its zip manifest. |
| `ensembles/ens_both.py` | Seed / architecture ensembling of stored model posteriors, for both modalities. Aligns members on clip path and repairs double-softmaxed posteriors before averaging. |
| `ensembles/ens_sel.py` | Selective ensembling — does curating the member pool help? |
| `ensembles/ens_err.py` | Splits error into the part members share and the part averaging can recover. |
| `ensembles/align_modalities.py` | Builds the clip intersection between two modalities' validation sets. |
| `preictal/` | Pre-ictal pipeline: lead-seizure cohort construction, hour-of-day-matched controls, horizon-wise effect sizes, circadian permutation null, post-ictal recovery. |
| `method/` | The methodology briefs the analyses were run against, including the known traps. |

## Findings

**The export is genuine.** Recomputing features from the source samples reproduces the
exported values to ~1e-5 relative — which is provably the 16-bit requantisation floor of
the local clip files, not error in the export (log-log fit of residual against half-LSB has
slope 1.0003, r = 0.966; clips whose sample grid happens to be lossless reproduce to <1e-12).

**Five bands carry about two numbers.** Within animal the bands correlate 0.26–0.89, one
component takes a median 71.4% of the variance, and the participation ratio is 1.83. A 1/f
generator with a level and a slope and *no band coupling* reproduces the entire correlation
structure. The pipeline itself contributes +0.0010 ± 0.0004 to adjacent-band correlation.

**The bands detect seizures but do not grade them.** Gamma rises 2.90 interictal SD during a
seizure and beta 2.23, while delta (+0.10) and theta (+0.44) have intervals crossing zero.
Across Racine stages 2→5 the gamma shift spans 0.23 SD. Within-session severe-vs-mild AUROC
is 0.46–0.55 for every band and for a 13-feature model.

**No pre-ictal change is detectable.** At every horizon from 60 min to 1 min before onset,
effects lie within ±0.24 SD, every animal-clustered interval spans zero, and nothing survives
Benjamini–Hochberg. The analysis is powered to ~0.35 SD, so this is a null with teeth for
group effects above that, and uninformative below ~0.3 SD.

**Ensembling helps video, not EEG.** Averaging stored video posteriors gains +0.003 / +0.011 /
+0.014 macro-F1 at detection / 3-class / 5-class for zero GPU. The same operation on EEG loses
to its own best single model under every honest selection rule.

## Running it

Python 3.11+ with `numpy` and `scipy`; `scikit-learn` only for the ensemble scripts.
Point the code at a data root:

```bash
export EEG_ROOT=/path/to/the/analysis/tree
python edf/provenance.py --animal RN244 --clips-per-day 1
python ensembles/ens_both.py
```

Every script defaults `EEG_ROOT` to the original analysis box, so it must be set elsewhere.

## What is not here

The data. The epoch export (~7 GB), the clip corpus (~36 GB) and the source EDFs are not in
this repository and are not redistributable from here. The exploratory scratch from the
adversarial verification passes is also omitted — each finding below was independently
re-measured by a second implementation, and only the corrected numbers are reported.
