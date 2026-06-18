# Nagauta Respiratory Coordination

This repository contains analysis code and data for:

Acting together beyond synchrony: Event-dependent alignment of breathing in an expert musical ensemble

---

## Dependencies

- MATLAB R2024b or later
- Signal Processing Toolbox (for `lowpass`, `downsample`)
- Statistics and Machine Learning Toolbox (for `ttest`, `corr`, `nchoosek`, `perms`)

No additional third-party toolboxes are required.

---

## Directory Structure

The scripts assume the following directory structure:

```
(root)/
├── Mfile/          # Analysis scripts (this repository)
└── Mat_file/
    ├── Hexoskin_data.mat
    ├── SectionTimepoint.mat
    └── MusicalEvent.mat
```

---

## Script Overview

Scripts are prefixed with `A` and can be run independently in any order.

### A — Respiratory Coordination Analyses

| Script | Description |
|--------|-------------|
| `A01_dtw_self_vs_other.m` | Compute normalised DTW distance between same-role (self-pairing) and different-role (other-pairing) respiration signals across takes; primary inference via exhaustive label-permutation test (6! = 720 permutations); paired t-test and Cohen's *d*z reported as descriptive measures |
| `A02_interperformer_correlation.m` | Compute lag-0 Pearson correlation (Fisher z) for all 15 performer dyads across sections and takes; primary inference via circular-shift surrogate test (1,000 iterations); one-sample t-test reported as descriptive measure |
| `A03_event_locked_respiration.m` | Extract and average thoracic respiration in a ±10 s window around each event type (Decel, Komi, Section); assess modulation against a 95% pointwise surrogate envelope from 1,000 random-event-time surrogates |

---

## Execution Order

All three scripts are independent and can be run in any order.

```
A01_dtw_self_vs_other
A02_interperformer_correlation
A03_event_locked_respiration
```

---

## Data

The following data files are included in `Mat_file/`. All signals are sampled at 128 Hz. All event and section times are in milliseconds relative to the start of each take.

| File | Variable | Description |
|------|----------|-------------|
| `Hexoskin_data.mat` | `Nagauta` | Thoracic respiration, ECG, and 3-axis accelerometer signals for 6 performers across 2 takes. Structure: `Nagauta.S1.{role}.{signal}` and `Nagauta.S2.{role}.{signal}`, where role is one of `Uta`, `Shamisen1`, `Kotsuzumi`, `Fue`, `Taiko`, `Shamisen2`, and signal is one of `respiration_thoracic`, `ECG`, `Acceleration`. |
| `SectionTimepoint.mat` | `Timepoint` | Section boundary times in ms; matrix of size [9 × 2] (rows: 8 section start times + 1 end time; columns: Take 1 / Take 2). Sections are labelled A–H. |
| `MusicalEvent.mat` | `Event` | Event onset times in ms for three event types. Structure: `Event.{type}.S1` and `Event.{type}.S2`, where type is one of `Decel` (tempo deceleration onsets), `Komi` (preparatory breath-cue events), or `Section` (section boundary onsets). |

---

## Notes

- All scripts use relative paths based on the `Mfile/` working directory. Before running, set the MATLAB working directory to `Mfile/`.
- The `A02` surrogate test uses a fixed random seed (`rng(20260205)`) for reproducibility.
- The 15 dyads in `A02` are formed from 6 performers and are therefore not statistically independent. The circular-shift surrogate test is the primary inferential test; the one-sample t-test is reported for descriptive purposes only (see manuscript for details).
- To analyse additional event types in `A03`, add the event type name to the `event_types` variable at the top of the script. The type must exist as a field in `MusicalEvent.mat`.

---

## Citation

If you use this code or data, please cite:

> (2026). Acting together beyond synchrony: Event-dependent alignment of breathing in an expert musical ensemble. *Acta Psychologica*. [DOI]

---

## License

This code is released under the MIT License.
