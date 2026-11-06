# Changelog

## Unreleased

- Support both photutils centroid column names (`x_centroid`/`y_centroid` and the older `xcentroid`/`ycentroid`), removing the deprecation warning with photutils >= 3.0.

## 0.1.0 (2026-11-05)

- Initial project structure: `src/banding_correction` package (`core.py`, `cli.py`), `banding-correction` command and `python -m banding_correction`.
- Algorithm taken unchanged from the original `banding_correction.py` script (output identical pixel for pixel on `bias_048.fits`).
- Synthetic tests (Bayer patterns, step detection, step size estimation, image correction).
- Apache-2.0 license.
