# Brief 2: the Mondino reference paper, and the physiology check it implies

Read `BRIEF.md` first for the export's mechanics, tools and gotchas. This file adds what is new.

## Answers the lab has now given (treat as authoritative)
* **Timezone:** "both time zones are the same" — the recording clock does not change across the 2023-11-05
  DST boundary. This matches what we measured (a single fixed-offset, strictly monotone clock,
  `epoch_start_native == recording_start + 10*epoch_index`, zero exceptions, 24 normal hours on 11-05).
  **Consequence: hour-of-day in `epoch_start_native` is internally consistent for the whole study and a
  circadian analysis on it is legitimate.** What is still unknown is the facility's lights-on / lights-off
  clock time, which is what anchors hour-of-day to the actual light schedule.
* **RN201 and RN203 have no seizures.** Their zero-annotation status is real, not a failed join. Their
  3,088.2 h are genuine seizure-free time and may be pooled as negatives.
* **The 3-channel `RoomC/DONE_RN203/11-15-2023.csv` export is intended.** channel_index 2 and 3 are
  deliberate. Analyses must still filter `channel_index == 1`; channel 3 is the 5 Hz high-passed copy.
* **`RoomD/RN245/` is to be discarded.** Drop those 5 files / 23 rows. This removes the export's only
  `channel_index == 0` rows and both ~1e9 uV2 amplitude outliers.
* **Spikes blanking is unexplained** — the lab does not know why every Spikes epoch has `features_valid=0`.
  It stays an open caveat.

## The reference paper
Mondino A, Cavelli M, Gonzalez J, Osorio L, Castro-Zaballa S, Costa A, Vanini G, Torterolo P.
"Power and Coherence in the EEG of the Rat: Impact of Behavioral States, Cortical Area, Lateralization and
Light/Dark Phases." *Clocks & Sleep* 2020;2(4):530-542. DOI 10.3390/clockssleep2040039. PMID 33317018.
PMCID PMC7768537.

**It is not an epilepsy paper.** Eleven normal Wistar males, no seizure model. Verbatim methods:
* Bands: "delta, 1-4 Hz; theta, 5-9 Hz; sigma, 10-15 Hz; beta, 16-30 Hz; low gamma (LG), 31-48 Hz;
  high gamma (HG) 52-95 Hz; high frequency oscillations (HFO), 105-200 Hz" — **non-contiguous**, with
  deliberate gaps at 4-5, 9-10, 48-52 and 95-105 Hz.
* Spectral method: "pwelch ... window = 30 s, noverlap = [], nfft = 2048, fs = 1024 ... 30 s sliding windows
  with half window overlap with a 0.5 Hz resolution."
* Relative power: "the absolute power of a specific frequency band/the sum of the power from 0.5 to 200 Hz."
* Recording: 8 intracranial screw electrodes (bilateral motor, somatosensory, visual cortex) + bipolar neck
  EMG, "amplified (x1000), filtered (0.1-500 Hz)", 1024 Hz 16 bit, Spike2.
* Housing: "12:12 h light/dark cycle (lights on at 6 a.m.)"; recordings "during the light (9 a.m. to 3 p.m.)
  and dark periods (9 p.m. to 3 a.m.)".
* States: W, light sleep + slow wave sleep grouped as NREM, and REM, scored from frontal/occipital EEG + EMG.
* Light/dark results: during **NREM**, "beta and LG power were larger during the dark than during the light
  phase in M1", and a cluster analysis found "sigma, beta, LG, HG and HFO bands were larger in the dark
  phase"; during **REM**, a theta cluster "higher during the day than during the night in V2";
  **no light/dark difference was detected during wakefulness**.
* State results: "delta, theta and sigma power during NREM sleep were significantly higher than during W and
  REM sleep"; "LG, HG and HFO powers were higher during W than during NREM in all the cortical areas";
  during REM "the relative weight of the theta band is highlighted".

## Why this matters to us, and what must be tested rather than assumed
Four hard non-comparabilities: different band edges, a different relative-power denominator (0.5-200 Hz vs our
0.5-80 Hz), a different window (30 s vs our 4 s inside a 10 s epoch), and a completely different electrode
type (intracranial screws over named cortical areas vs a single DSI telemetry biopotential channel).
Absolute power magnitudes are therefore **not** comparable; only directions of effect might transfer.

The load-bearing warning: the paper's light/dark spectral effects live almost entirely **inside NREM**, with
**none in wakefulness**. Our export has **no sleep staging at all**. So an unstaged light/dark band-power
contrast on our data mostly measures the light/dark difference in how much the animal *sleeps*, not a
within-state spectral change. This project has already been burned once by exactly this class of error - the
EEG severity AUROC that turned out to be session base rate (see the settled finding that EEG carries
detection but not within-event severity). Do not repeat it.

## Report format
Same structured format as Brief 1. Quote exact numbers, say what was exhaustive vs sampled, and flag any
claim that a state confound could explain.
