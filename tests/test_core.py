import numpy as np
import pytest

from banding_correction.core import (bayer_slices, compute_row_offsets,
                                     correct_image, estimate_step)


@pytest.mark.parametrize("pattern,expected", [
    ("RGGB", {"R": (0, 0), "G1": (0, 1), "G2": (1, 0), "B": (1, 1)}),
    ("BGGR", {"B": (0, 0), "G1": (0, 1), "G2": (1, 0), "R": (1, 1)}),
    ("GRBG", {"G1": (0, 0), "R": (0, 1), "B": (1, 0), "G2": (1, 1)}),
    ("GBRG", {"G1": (0, 0), "B": (0, 1), "R": (1, 0), "G2": (1, 1)}),
])
def test_bayer_slices(pattern, expected):
    assert bayer_slices(pattern) == expected


def test_bayer_slices_invalid():
    with pytest.raises(ValueError):
        bayer_slices("XXXX")


def synthetic_profile(n=2000, step=1.3, sigma=0.07, seed=0):
    rng = np.random.default_rng(seed)
    true = np.zeros(n)
    rows = np.arange(100, n - 50, 150)
    level = 0.0
    for i, r in enumerate(rows):
        level += step if i % 2 == 0 else -step
        true[r:] = level
    return true + rng.normal(0, sigma, n), true


def test_row_offsets_flatten_profile():
    profile, true = synthetic_profile()
    offsets, _, _ = compute_row_offsets(profile, 0.07, step=1.3)
    residual = profile + offsets
    assert np.nanstd(residual[16:]) < 0.15


def test_estimate_step():
    profile, _ = synthetic_profile(step=1.3)
    assert estimate_step(profile, 0.07) == pytest.approx(1.3, abs=0.1)


def test_estimate_step_flat_channel():
    rng = np.random.default_rng(1)
    assert estimate_step(rng.normal(0, 0.07, 2000), 0.07) is None


def test_correct_image_removes_banding():
    rng = np.random.default_rng(2)
    h, w = 800, 3000  # full image; each channel is 400 x 1500
    base = rng.normal(1000, 4, (h, w)).astype(np.float32)
    banded = base.copy()
    level = np.zeros(h // 2)
    for r in range(40, h // 2 - 20, 60):
        level[r:] += 1.3 if (r // 60) % 2 == 0 else -1.3
    banded += np.repeat(level, 2)[:, None].astype(np.float32)
    corrected, _, steps = correct_image(banded, "RGGB", verbose=False)
    err_before = np.abs(banded.mean(axis=1) - base.mean(axis=1)).mean()
    err_after = np.abs(corrected.mean(axis=1) - base.mean(axis=1)).mean()
    assert err_after < err_before / 3
