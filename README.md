# Nagauta Respiratory Synchrony

This repository contains analysis code and data for:

Acting together beyond synchrony: Event-dependent alignment of breathing in an expert musical ensemble

---

## Dependencies

- MATLAB R2024b or later
- Signal Processing Toolbox (for `lowpass`, `downsample`)
- Statistics and Machine Learning Toolbox (for `ttest`, `corr`, `nchoosek`)

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

### A — Respiratory Synchrony Analyses

| Script | Description |
|--------|-------------|
| `A01_dtw_self_vs_other.m` | Compute normalised DTW distance between same-role (self-pairing) and different-role (other-pairing) respiration signals across takes; compare with a paired t-test |
| `A02_interperformer_correlation.m` | Compute lag-0 Pearson correlation (Fisher z) for all 15 performer dyads across sections and takes; test against chance using a circular-shift surrogate test (primary) and a one-sample t-test (descriptive) |
| `A03_event_locked_respiration.m` | Extract and average thoracic respiration in a ±10 s window around each event type (Decel, Komi, Section); compute a 95% pointwise surrogate envelope from 1000 random-event-time surrogates |

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

The following data files are included in `Mat_file/`:

| File | Variable | Description |
|------|----------|-------------|
| `Hexoskin_data.mat` | `Nagauta` | Thoracic respiration, ECG, and accelerometer signals (128 Hz) for 6 performers across 2 takes |
| `SectionTimepoint.mat` | `Timepoint` | Section boundary times in ms; matrix of size [9 × 2] (8 sections + end point, 2 takes) |
| `MusicalEvent.mat` | `Event` | Event onset times in ms for Decel, Komi, and Section events, stored separately for each take (S1, S2) |

See `DATA_DESCRIPTION.md` for the full data structure.

---

## Notes

- All scripts use relative paths based on the `Mfile/` working directory. Before running, set the MATLAB working directory to `Mfile/`.
- The `A02` surrogate test uses a fixed random seed (`rng(20260205)`) for reproducibility.
- The 15 dyads in `A02` are formed from 6 performers and are therefore not statistically independent. The circular-shift surrogate test is the primary inferential test; the one-sample t-test is reported for descriptive purposes only (see manuscript for details).
- To analyse additional event types in `A03`, add the event type name to the `event_types` variable at the top of the script. The type must exist as a field in `MusicalEvent.mat`.

---

## Citation

If you use this code or data, please cite:

> ([Year]). Acting together beyond synchrony: Event-dependent alignment of breathing in an expert musical ensemble. *Acta Psychologica*. [DOI]

---

## License

This code is released under the MIT License.
