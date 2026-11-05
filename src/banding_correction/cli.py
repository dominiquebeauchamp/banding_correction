"""Command line interface of banding_correction."""
import argparse
import os
import sys

from astropy.nddata import CCDData

from .core import bayer_slices, correct_image, plot_profiles


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Remove horizontal banding steps from a CFA (Bayer) FITS image.')
    p.add_argument('input', help='FITS image to correct')
    p.add_argument('output', help='corrected FITS image to write')
    p.add_argument('--fwhm', type=float, default=3.0,
                   help='star FWHM in channel pixels (default: %(default)s)')
    p.add_argument('--detect-sigma', type=float, default=5.0,
                   help='star detection threshold in units of background std (default: %(default)s)')
    p.add_argument('--mask-radius-factor', type=float, default=3.0,
                   help='star mask radius = factor * FWHM (default: %(default)s)')
    p.add_argument('--clip-sigma', type=float, default=3.0,
                   help='sigma clipping of each row when averaging (default: %(default)s)')
    p.add_argument('--window', type=int, default=16,
                   help='number of previous rows averaged to detect a step (default: %(default)s)')
    p.add_argument('--confirm-rows', type=int, default=3,
                   help='number of rows averaged ahead of the current one to confirm a step '
                        '(default: %(default)s)')
    p.add_argument('--threshold-sigma', type=float, default=4.5,
                   help='step detection threshold, in standard deviations of the difference '
                        'between the confirmation rows and the previous window (default: %(default)s)')
    p.add_argument('--step', type=float, default=None,
                   help='size in ADU of the row steps; by default it is estimated '
                        'automatically for each channel')
    p.add_argument('--default-step', type=float, default=1.0,
                   help='step used when no channel has detectable steps and --step is not '
                        'given (default: %(default)s)')
    p.add_argument('--overwrite', action='store_true',
                   help='overwrite the output file if it exists')
    p.add_argument('--plot', metavar='PNG',
                   help='save a plot of original/corrected vertical profiles to this file')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if not os.path.isfile(args.input):
        sys.exit(f"error: input file not found: {args.input}")
    if os.path.exists(args.output) and not args.overwrite:
        sys.exit(f"error: {args.output} already exists (use --overwrite)")
    if args.window < 2:
        sys.exit("error: --window must be >= 2")
    if args.confirm_rows < 1:
        sys.exit("error: --confirm-rows must be >= 1")

    ccd = CCDData.read(args.input, unit='adu')
    if ccd.data.ndim != 2:
        sys.exit(f"error: expected a 2D image, got shape {ccd.data.shape}")

    pattern = str(ccd.header.get('BAYERPAT', 'RGGB')).strip().upper()
    try:
        bayer_slices(pattern)
    except ValueError as exc:
        sys.exit(f"error: {exc}")
    print(f"Bayer pattern: {pattern}")

    corrected, profiles, steps = correct_image(
        ccd.data, pattern, fwhm=args.fwhm, detect_sigma=args.detect_sigma,
        radius_factor=args.mask_radius_factor, clip_sigma=args.clip_sigma,
        window=args.window, confirm_rows=args.confirm_rows,
        threshold_sigma=args.threshold_sigma, step=args.step,
        default_step=args.default_step)

    out = ccd.copy()
    out.data = corrected
    out.header['HISTORY'] = 'Banding correction (banding_correction package)'
    out.header['HISTORY'] = (f'window={args.window} confirm_rows={args.confirm_rows} '
                             f'threshold_sigma={args.threshold_sigma} clip_sigma={args.clip_sigma}')
    out.header['HISTORY'] = 'steps (ADU): ' + ' '.join(f'{k}={v:.2f}' for k, v in steps.items())
    out.header['HISTORY'] = (f'stars: fwhm={args.fwhm} detect_sigma={args.detect_sigma} '
                             f'radius_factor={args.mask_radius_factor}')
    out.write(args.output, overwrite=args.overwrite)
    print(f"Wrote {args.output}")

    if args.plot:
        plot_profiles(profiles, args.plot)
        print(f"Wrote {args.plot}")


if __name__ == '__main__':
    main()
