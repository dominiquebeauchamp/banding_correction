"""Remove horizontal banding (row-wise ~1 ADU steps) from a CFA (Bayer) FITS image.

Algorithm (see banding_002.ipynb and banding_003.ipynb):
  1. Split the image into its four Bayer channels (layout given by BAYERPAT).
  2. Detect stars on the mean of the two green channels and mask them (NaN) in
     all channels, using the same mask.
  3. For each channel, compute the vertical profile (clipped mean of each row).
  4. Estimate the profile noise from the MAD of its first differences.
  5. Walk down the profile, keeping a cumulative offset. At each row, the mean of
     the next `confirm_rows` values (plus the current offset) is compared with the
     mean of the previous `window` corrected values. A deviation larger than
     `threshold_sigma` standard deviations of that difference signals a step.
     Averaging over several rows makes the test robust to isolated noise spikes.
  5b. The step position is then refined: among the rows covered by the test, the
     one that best splits the profile into "no step / step of `step` ADU" is
     chosen, and the offset is applied from that row on.
  5c. The step size is estimated automatically for each channel (unless --step is
     given): it is the value for which the corrected profile is the flattest.
  6. Add the per-row offset to every pixel of the corresponding rows of the
     original (unmasked) channel and write the result.
"""
import warnings

import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from scipy.ndimage import distance_transform_edt

VALID_PATTERNS = ('RGGB', 'BGGR', 'GRBG', 'GBRG')


# --------------------------------------------------------------------------- #
# Bayer helpers
# --------------------------------------------------------------------------- #
def bayer_slices(pattern):
    """Return {channel_name: (row_offset, col_offset)} for a Bayer pattern.

    The pattern string lists the colours at positions (0,0), (0,1), (1,0), (1,1)
    of the array. The two greens are named G1 and G2 in reading order.
    """
    pattern = pattern.strip().upper()
    if pattern not in VALID_PATTERNS:
        raise ValueError(f"Unsupported BAYERPAT '{pattern}' "
                         f"(expected one of {', '.join(VALID_PATTERNS)})")
    positions = [(0, 0), (0, 1), (1, 0), (1, 1)]
    slices = {}
    n_green = 0
    for letter, pos in zip(pattern, positions):
        if letter == 'G':
            n_green += 1
            slices[f'G{n_green}'] = pos
        else:
            slices[letter] = pos
    return slices


def extract_channel(image, offsets):
    """Sub-sample one Bayer channel from the full image."""
    dy, dx = offsets
    return image[dy::2, dx::2]


# --------------------------------------------------------------------------- #
# Star masking
# --------------------------------------------------------------------------- #
def detect_star_mask(green_mean, fwhm, detect_sigma, radius_factor):
    """Return a boolean mask (True = star) on the channel grid.

    Stars are detected with DAOStarFinder on the mean of the green channels;
    each centroid is surrounded by a disc of radius `radius_factor * fwhm`.
    """
    _, median, std = sigma_clipped_stats(green_mean, sigma=3.0)
    finder = DAOStarFinder(fwhm=fwhm, threshold=detect_sigma * std)
    sources = finder(green_mean - median)
    if sources is None or len(sources) == 0:
        return np.zeros(green_mean.shape, dtype=bool), 0

    centers = np.zeros(green_mean.shape, dtype=bool)
    x = np.clip(np.round(sources['xcentroid']).astype(int), 0, green_mean.shape[1] - 1)
    y = np.clip(np.round(sources['ycentroid']).astype(int), 0, green_mean.shape[0] - 1)
    centers[y, x] = True
    mask = distance_transform_edt(~centers) <= radius_factor * fwhm
    return mask, len(sources)


# --------------------------------------------------------------------------- #
# Vertical profile and step correction
# --------------------------------------------------------------------------- #
def vertical_profile(channel, sigma=3.0):
    """Clipped mean of each row (NaN ignored).

    The row median is only used to centre the clipping. The profile is a mean
    rather than a median to avoid the quantisation of integer ADU data.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)  # all-NaN rows
        ref = np.nanmedian(channel, axis=1, keepdims=True)
        s = 1.4826 * np.nanmedian(np.abs(channel - ref), axis=1, keepdims=True)
        clipped = np.where(np.abs(channel - ref) <= sigma * s, channel, np.nan)
        return np.nanmean(clipped, axis=1)


def estimate_profile_noise(profile):
    """Robust noise estimate of the profile from its first differences.

    Steps are rare outliers in the differences, so the MAD ignores them. The
    difference of two independent values has twice the variance, hence 1/sqrt(2).
    """
    dp = np.diff(profile)
    mad = np.nanmedian(np.abs(dp - np.nanmedian(dp)))
    return 1.4826 * mad / np.sqrt(2)


def compute_row_offsets(profile, sigma_p, window=16, confirm_rows=3,
                        threshold_sigma=4.5, step=1.0):
    """Cumulative-offset step correction of a vertical profile.

    Starting at row `window`, the mean of the `confirm_rows` rows y..y+m-1 (shifted
    by the current offset) is compared with the mean of the `window` previous
    corrected rows. The difference of two means has a standard deviation of
    sigma_p * sqrt(1/m + 1/window); a deviation beyond `threshold_sigma` times that
    value signals a step. Compared with a single-row test, a real step of 1 ADU is
    seen with a much larger significance, so the threshold can be raised and isolated
    noise spikes no longer create false steps.

    The detection can fire up to m-1 rows before the true step, because the test
    already sees a few rows on the other side. The position is therefore refined:
    among the candidate rows y..y+m-1, the one minimising the squared error between
    the rows y..y+2m-1 (relative to the reference mean) and the model "0 before the
    step, +-step after it" is chosen, and the offset is applied from that row on.
    A rise of the profile is compensated by subtracting `step`, a fall by adding it.
    NaN rows are skipped but keep the current offset.

    Returns (offsets, n_down, n_up): one offset per row, and the number of
    corrections of -step and +step.
    """
    m = max(1, int(confirm_rows))
    loc = 2 * m                     # number of rows used to locate the step
    n = len(profile)
    corrected = profile.copy()
    offsets = np.zeros(n)
    offset = 0.0
    pending = None                  # (row, direction): step located, not yet applied
    n_down = n_up = 0
    for y in range(window, n):
        if pending is not None and y == pending[0]:
            offset -= pending[1] * step
            if pending[1] > 0:
                n_down += 1
            else:
                n_up += 1
            pending = None
        if not np.isnan(profile[y]) and pending is None:
            previous = corrected[y - window:y]
            ahead = profile[y:y + m] + offset
            n_prev, n_ahead = np.sum(~np.isnan(previous)), np.sum(~np.isnan(ahead))
            if n_prev >= 2 and n_ahead >= 1:
                ref_mean = np.nanmean(previous)
                delta = np.nanmean(ahead) - ref_mean
                threshold = threshold_sigma * sigma_p * np.sqrt(1.0 / n_ahead + 1.0 / n_prev)
                if abs(delta) > threshold:
                    direction = 1.0 if delta > 0 else -1.0
                    # Locate the step among the rows covered by the test
                    u = profile[y:y + loc] + offset - ref_mean
                    idx = np.arange(len(u))
                    best_s, best_cost = 0, np.inf
                    for s in range(min(m, len(u))):
                        model = np.where(idx >= s, direction * step, 0.0)
                        cost = np.nansum((u - model) ** 2)
                        if cost < best_cost:
                            best_s, best_cost = s, cost
                    if best_s == 0:
                        offset -= direction * step
                        if direction > 0:
                            n_down += 1
                        else:
                            n_up += 1
                    else:
                        pending = (y + best_s, direction)
        if not np.isnan(profile[y]):
            corrected[y] = profile[y] + offset
        offsets[y] = offset
    return offsets, n_down, n_up


def _profile_spread(profile, sigma_p, step, window, confirm_rows, threshold_sigma):
    """Robust spread (MAD) of the corrected profile, and number of offset changes."""
    offsets, n_down, n_up = compute_row_offsets(
        profile, sigma_p, window=window, confirm_rows=confirm_rows,
        threshold_sigma=threshold_sigma, step=step)
    corrected = profile + offsets
    spread = 1.4826 * np.nanmedian(np.abs(corrected - np.nanmedian(corrected)))
    return spread, n_down + n_up


def estimate_step(profile, sigma_p, window=16, confirm_rows=3, threshold_sigma=4.5,
                  step_range=(0.6, 2.5), min_events=4):
    """Estimate the size (ADU) of the row steps of a vertical profile.

    The row offsets are not necessarily integer: they are applied before the ADC,
    and the ~thousands of pixels averaged in each row dither the quantisation, so the
    profile can jump by any amount (e.g. 1.3 ADU). If the correction step is wrong,
    part of each step remains uncorrected and the residual accumulates.

    The correction is run for a range of candidate steps (coarse grid, then a finer
    grid around the best one) and the step giving the flattest corrected profile
    (smallest robust spread) is kept. A channel without steps has a flat curve: if
    fewer than `min_events` offset changes are found with a step of 1 ADU, the channel
    is considered step-free and None is returned.
    """
    def run(st):
        return _profile_spread(profile, sigma_p, st, window, confirm_rows, threshold_sigma)

    if run(1.0)[1] < min_events:
        return None
    coarse = np.arange(step_range[0], step_range[1] + 1e-9, 0.15)
    best = coarse[int(np.argmin([run(st)[0] for st in coarse]))]
    fine = np.arange(max(best - 0.15, 0.1), best + 0.15 + 1e-9, 0.03)
    scores = [run(st)[0] for st in fine]
    # Mean of the candidates within 2% of the minimum: less sensitive to noise
    # than a bare argmin on a flat-bottomed curve.
    good = np.array(scores) <= min(scores) * 1.02
    return float(np.mean(fine[good]))


# --------------------------------------------------------------------------- #
# Whole-image correction
# --------------------------------------------------------------------------- #
def correct_image(data, pattern, fwhm=3.0, detect_sigma=5.0, radius_factor=3.0,
                  clip_sigma=3.0, window=16, confirm_rows=3, threshold_sigma=4.5,
                  step=None, default_step=1.0, verbose=True):
    """Correct row banding in a CFA image.

    `step` is the size in ADU of the row steps. If None, it is estimated for each
    channel with `estimate_step`; a channel without detectable steps uses the median
    of the other channels' estimates, or `default_step` if no channel has any.

    Returns (corrected_image, profiles, steps): profiles maps each channel name to
    (original_profile, corrected_profile) for optional plotting, and steps maps each
    channel name to the step size used.
    """
    slices = bayer_slices(pattern)
    data = data.astype(np.float32)
    corrected = data.copy()

    # Analysis grid: crop to a common shape in case the image has odd dimensions.
    h, w = data.shape[0] // 2, data.shape[1] // 2
    channels = {name: extract_channel(data, off)[:h, :w] for name, off in slices.items()}

    # Star detection on the mean of the green channels, shared by all channels.
    green_mean = (channels['G1'] + channels['G2']) / 2.0
    star_mask, n_stars = detect_star_mask(green_mean, fwhm, detect_sigma, radius_factor)
    if verbose:
        print(f"{n_stars} stars detected, {star_mask.mean() * 100:.1f}% of pixels masked")

    # Vertical profile and noise of each channel (stars masked)
    raw_profiles, sigmas = {}, {}
    for name in slices:
        masked = channels[name].copy()
        masked[star_mask] = np.nan
        raw_profiles[name] = vertical_profile(masked, sigma=clip_sigma)
        sigmas[name] = estimate_profile_noise(raw_profiles[name])

    # Step size of each channel
    steps, origin = {}, {}
    if step is not None:
        steps = {name: float(step) for name in slices}
        origin = {name: 'given' for name in slices}
    else:
        estimates = {name: estimate_step(raw_profiles[name], sigmas[name], window=window,
                                         confirm_rows=confirm_rows,
                                         threshold_sigma=threshold_sigma)
                     for name in slices}
        found = [v for v in estimates.values() if v is not None]
        fallback = float(np.median(found)) if found else default_step
        for name, value in estimates.items():
            steps[name] = value if value is not None else fallback
            origin[name] = 'estimated' if value is not None else (
                'from other channels' if found else 'default')

    profiles = {}
    for name, off in slices.items():
        profile, sigma_p, st = raw_profiles[name], sigmas[name], steps[name]
        offsets, n_down, n_up = compute_row_offsets(
            profile, sigma_p, window=window, confirm_rows=confirm_rows,
            threshold_sigma=threshold_sigma, step=st)
        if verbose:
            print(f"{name}: sigma_p = {sigma_p:.4f} ADU, step = {st:.2f} ADU ({origin[name]}), "
                  f"{n_down} x -step, {n_up} x +step, final offset = {offsets[-1]:+.2f} ADU")

        # Apply the offset to the whole rows of the original channel (stars included).
        # Odd image sizes leave one extra row in some channels: reuse the last offset.
        target = corrected[off[0]::2, off[1]::2]
        full_rows = np.pad(offsets, (0, target.shape[0] - len(offsets)), mode='edge')
        target += full_rows[:, None].astype(np.float32)  # view: modifies `corrected`

        profiles[name] = (profile, profile + offsets)

    return corrected, profiles, steps


def plot_profiles(profiles, path):
    """Save original/corrected vertical profiles of each channel to an image."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(profiles), 1, figsize=(10, 3 * len(profiles)), sharex=True)
    for ax, (name, (orig, corr)) in zip(np.atleast_1d(axes), profiles.items()):
        ax.plot(orig, 'b-', linewidth=0.8, label='original')
        ax.plot(corr, 'r-', linewidth=0.6, label='corrected')
        ax.set_ylabel(f'{name} (ADU)')
        ax.grid()
        ax.legend()
    np.atleast_1d(axes)[-1].set_xlabel('row')
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
