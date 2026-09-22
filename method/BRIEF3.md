# Brief 3: the pre-ictal question

Read `BRIEF.md` then `BRIEF2.md` first. This adds the pre-ictal task and the traps specific to it.

## The question
Do the exported band features change **before** a seizure onset, and if so is that change usable to
*predict* one? These are two different claims and must be reported separately:
* **Characterisation** — do pre-ictal epochs differ from properly matched interictal epochs?
* **Prediction** — prospectively, out of sample, can an alarm be raised with a useful sensitivity at a
  tolerable false-alarm rate? This is the claim that matters and the one that is almost always overstated.

## Established, do not re-derive (challenge only with evidence)
* Base set: 11,722,320 epochs, 23 animals. Filters: `channel_index == 1`, `features_valid == 1`,
  exclude `RoomD/RN245/` and the `.interrupted` file. 10-second epochs.
* 13,356 annotation events. Resolve each id's stage ONLY from epochs where it is the sole id.
  Onset = start of the earliest epoch with `annotation_overlap_fraction > 0`.
* **Ictal effect sizes** (median z vs interictal, animal-clustered, n=19): gamma +2.90, beta +2.23,
  alpha +1.09, theta +0.44 (n.s.), delta +0.10 (n.s.), line_length +3.28. Detection AUROC within
  session: line_length 0.868, gamma 0.831, beta 0.781, alpha 0.704, theta 0.579, delta 0.546.
* **The bands are ~2 dimensions**, not 5: PC1 takes 71.4% of within-animal variance, participation
  ratio 1.83, and a 1/f level+slope generator reproduces the whole correlation structure. So test
  **level (log total power) and slope/tilt** as the primary axes and treat the five bands as secondary
  — it cuts the multiplicity problem honestly rather than testing five correlated things.
* **Seizures cluster in the LIGHT/rest phase**: 55.9 onsets per 100 recorded h in 06:00-18:00 vs 34.5
  at night, ratio 1.62 (animal-clustered CI 1.50-1.80), 20/21 animals.
* **Band power has a strong light-locked rhythm**: alpha swings 0.64-1.30 of its own daily median,
  stepping at 06:00 and ~17:50. line_length runs the OTHER way (higher at night).
* `features_valid = 0` on 2.53% of rows, and the blanking is **diurnal** (dark 3.34% vs light 1.73%)
  and **severity-correlated**. It is a selection channel, not noise.
* **Spikes events carry no features at all** (100% blanked) and cannot appear.
* **RN201 and RN203 are confirmed by the lab to have NO seizures** — 3,088 h of clean seizure-free
  recording. They are the correct place to measure a false-alarm rate with no censoring.

## The five traps. Every one of these has produced a published false positive somewhere.
1. **Circadian leakage — the biggest one.** Seizures are 1.62x more likely in the light phase and band
   power swings ~2x with the light cycle. A pre-ictal window is therefore *more likely to be daytime*
   than a randomly drawn interictal control, so an uncontrolled comparison will find a "pre-ictal
   signature" that is pure time-of-day. **Controls must be matched on hour-of-day within animal**, and
   you must also report the result after regressing out the hour-of-day profile.
2. **Seizure clustering — post-ictal masquerading as pre-ictal.** Rodent seizures cluster. The window
   before seizure N is often the *post-ictal* period of seizure N-1. Restrict the primary analysis to
   **lead seizures** (no other event within at least the horizon + a margin, state your choice) and
   report how many events survive that filter. Measure the post-ictal recovery time course explicitly
   so the reader can see how far contamination reaches.
3. **Onset-time imprecision.** The annotated onset is not exact and early ictal activity may precede
   it. Report every horizon separately and treat the last minute before onset as suspect, not as
   evidence of prediction.
4. **Selection by blanking.** `features_valid = 0` is diurnal and severity-correlated, so dropping
   those epochs is not missing-at-random. Report how many pre-ictal and control epochs are lost.
5. **Unit of inference.** Animal, never epoch. 11.7 M epochs will make anything "significant".
   Use leave-one-animal-out for anything predictive, and animal-clustered intervals throughout.

## Reporting rules
State the estimand before the number. Give effect sizes with animal-clustered intervals, not p-values
alone. Correct for multiplicity across bands x horizons and say how. If the honest answer is "no
detectable pre-ictal change", say that plainly — a clean null here is a genuinely useful result for
this project and is *far* more likely a priori than a positive. Do not report a prediction result
without a false-alarm rate and a comparison against a predictor that knows only the time of day.
