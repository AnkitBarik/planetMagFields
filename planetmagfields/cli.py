#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI tool for plotting planetary magnetic fields.

This module provides a command-line interface for visualizing magnetic field
data from planets in the solar system at various radial levels and with
customizable plot options.
"""

import sys
import json
import logging
import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from .libbfield import plotAllFields
from .planet import Planet
from . import __version__

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s"
)
logger = logging.getLogger(__name__)

DEFAULT_CONFIG_PATH = Path(__file__).parent / "config.json"

# Keys that are valid in the config file and their expected types
_CONFIG_TYPES = {
    "planet": str,
    "model": str,
    "r": float,
    "unit": str,
    "cmap": str,
    "levels": int,
    "proj": str,
    "vmin": float,
    "vmax": float,
    "save_path": str,
    "dark": bool,
    "verbose": bool,
    "trajectory": str,
}


def load_config(path=None):
    """Load options from a JSON config file.

    Checks, in order:
      1. The path passed via ``--config`` on the command line.
      2. The ``config.json`` bundled with the package (always present).

    Returns a dict of validated key/value pairs (unknown keys are ignored).
    """
    if path is None:
        if not DEFAULT_CONFIG_PATH.exists():
            return {}
        path = DEFAULT_CONFIG_PATH

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with path.open() as fh:
        try:
            raw = json.load(fh)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid JSON in config file {path}: {exc}") from exc

    if not isinstance(raw, dict):
        raise ValueError(f"Config file must contain a JSON object, got {type(raw).__name__}")

    config = {}
    for key, value in raw.items():
        if key.startswith("_"):  # reserved for comments / metadata
            continue
        if key not in _CONFIG_TYPES:
            logger.warning(f"Unknown config key '{key}' — ignored")
            continue
        expected = _CONFIG_TYPES[key]
        if value is None:
            config[key] = None
            continue
        try:
            config[key] = expected(value)
        except (TypeError, ValueError):
            raise ValueError(
                f"Config key '{key}' expects {expected.__name__}, got {type(value).__name__}"
            )

    return config


def setup_dark_background():
    """Configure matplotlib styling for dark background."""
    plt.style.use("dark_background")
    plt.rcParams.update({
        "axes.facecolor": "#1b1b1b",
        "figure.facecolor": "#1b1b1b",
        "figure.edgecolor": "#1b1b1b",
        "savefig.facecolor": "#1b1b1b",
        "savefig.edgecolor": "#1b1b1b"
    })


def validate_arguments(args):
    """Validate command-line arguments."""
    if args.r <= 0:
        raise ValueError("Radius must be positive")

    if args.levels < 1:
        raise ValueError("Number of levels must be at least 1")

    if args.vmin is not None and args.vmax is not None:
        if args.vmin > args.vmax:
            raise ValueError("vmin must be less than or equal to vmax")


def run_trajectory(args):
    """Compute field components along a trajectory defined by a CSV file."""
    import io
    import numpy as np

    traj_path = Path(args.trajectory)
    if not traj_path.exists():
        raise FileNotFoundError(f"Trajectory file not found: {traj_path}")

    # Detect whether the first row is a header by trying to parse it as floats.
    with traj_path.open() as fh:
        first_line = fh.readline().strip()

    has_header = False
    try:
        [float(v) for v in first_line.split(',')]
    except ValueError:
        has_header = True

    try:
        if has_header:
            data = np.genfromtxt(traj_path, delimiter=',', names=True, dtype=float)
            required = {'r', 'theta', 'phi'}
            missing = required - set(data.dtype.names)
            if missing:
                raise ValueError(
                    f"Trajectory CSV header is missing column(s): "
                    f"{', '.join(sorted(missing))}. Expected: r,theta,phi"
                )
            r, theta, phi = data['r'], data['theta'], data['phi']
        else:
            data = np.genfromtxt(traj_path, delimiter=',', dtype=float)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            if data.shape[1] < 3:
                raise ValueError(
                    "Trajectory CSV must have at least 3 columns (r, theta, phi)"
                )
            r, theta, phi = data[:, 0], data[:, 1], data[:, 2]
    except ValueError:
        raise
    except Exception as exc:
        raise ValueError(f"Could not read trajectory file: {exc}") from exc

    logger.info(f"Running orbit_path for {args.planet} over {len(r)} points...")
    planet = Planet(name=args.planet, model=args.model)
    planet.orbit_path(r, theta, phi)

    out = np.column_stack([r, theta, phi,
                           planet.br_orb, planet.btheta_orb, planet.bphi_orb])
    header = f"r,theta,phi,Br,Btheta,Bphi  (field unit: {planet.unit})"

    if args.save_path:
        out_path = Path(args.save_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(out_path, out, delimiter=',', header=header, comments='')
        logger.info(f"Results saved to {out_path}")
    else:
        buf = io.StringIO()
        np.savetxt(buf, out, delimiter=',', header=header, comments='')
        print(buf.getvalue(), end='')


def plot_all_planets(args):
    """Plot magnetic fields for all planets."""
    import importlib.resources as pkg_resources

    logger.info("Plotting magnetic fields for all planets...")

    data_dir = Path(pkg_resources.files("planetmagfields") / "data")

    plotAllFields(
        datDir=str(data_dir) + "/",
        r=args.r,
        levels=args.levels,
        cmap=args.cmap,
        proj=args.proj,
        unit=args.unit,
        vmin=args.vmin,
        vmax=args.vmax
    )

    plt.tight_layout()
    plt.subplots_adjust(
        top=0.895,
        bottom=0.035,
        left=0.023,
        right=0.976,
        hspace=0.38,
        wspace=0.109
    )


def plot_single_planet(args):
    """Plot magnetic field for a single planet."""
    logger.info(f"Plotting magnetic field for {args.planet}...")

    try:
        planet = Planet(name=args.planet, r=args.r, model=args.model)
        planet.plot(r=args.r, levels=args.levels, cmap=args.cmap, proj=args.proj)
    except Exception as e:
        logger.error(f"Failed to plot {args.planet}: {e}")
        raise


def save_plot(output_path):
    """Save the current plot to a file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logger.info(f"Plot saved to {output_path}")


def create_parser():
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        prog='magfield',
        description='Plot planetary magnetic fields from the solar system.',
        epilog='For more information, visit: https://ankitbarik.github.io/planetMagFields/'
    )

    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__version__}'
    )

    # Planet selection
    planet_group = parser.add_argument_group('Planet Selection')
    planet_group.add_argument(
        '-p', '--planet',
        type=str,
        default='earth',
        dest='planet',
        help='Planet name (default: earth)'
    )
    planet_group.add_argument(
        '-o', '--model',
        type=str,
        default=None,
        dest='model',
        help='Specific model to use (uses latest model by default)'
    )

    # Radial and data options
    data_group = parser.add_argument_group('Data Options')
    data_group.add_argument(
        '-r', '--radius',
        type=float,
        default=1.0,
        dest='r',
        help='Radial level scaled to planetary radius (default: 1.0)'
    )
    data_group.add_argument(
        '-u', '--unit',
        type=str,
        default='muT',
        dest='unit',
        help='Unit of magnetic field (default: muT)'
    )
    data_group.add_argument(
        '--trajectory',
        type=str,
        default=None,
        dest='trajectory',
        metavar='FILE',
        help=(
            'CSV file with columns r, theta, phi (co-latitude and longitude in radians). '
            'Computes Br, Btheta, Bphi along the path instead of plotting. '
            'Output is written to stdout or to --save FILE as a CSV.'
        )
    )

    # Plotting options
    plot_group = parser.add_argument_group('Plot Options')

    plot_group.add_argument(
        '--dark',
        action='store_true',
        dest='dark',
        help='Use dark background for plots'
    )
    plot_group.add_argument(
        '-c', '--cmap',
        type=str,
        default='RdBu_r',
        dest='cmap',
        help='Colormap to use (default: RdBu_r)'
    )
    plot_group.add_argument(
        '-l', '--levels',
        type=int,
        default=20,
        dest='levels',
        help='Number of contour levels (default: 20)'
    )
    plot_group.add_argument(
        '-m', '--mapproj',
        type=str,
        default='Mollweide',
        dest='proj',
        help='Map projection type (default: Mollweide)'
    )
    plot_group.add_argument(
        '--vmin',
        type=float,
        default=None,
        dest='vmin',
        help='Minimum value for colorscale'
    )
    plot_group.add_argument(
        '--vmax',
        type=float,
        default=None,
        dest='vmax',
        help='Maximum value for colorscale'
    )

    # Output options
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument(
        '-s', '--save',
        type=str,
        default=None,
        dest='save_path',
        help='Save plot to file instead of displaying (e.g., output.png)'
    )
    output_group.add_argument(
        '-v', '--verbose',
        action='store_true',
        dest='verbose',
        help='Enable verbose output'
    )
    output_group.add_argument(
        '--config',
        type=str,
        default=None,
        dest='config',
        metavar='FILE',
        help=(
            'Path to a JSON config file. '
            f'Defaults to the config.json bundled with the package. '
            'CLI flags override config values.'
        )
    )

    return parser


def main(argv=None):
    """Main entry point for the CLI."""
    parser = create_parser()

    # --- Config file: parse only --config / --verbose first so we can load
    #     the file before the full parse, then set its values as defaults so
    #     explicit CLI flags still win.
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument('--config', default=None)
    pre_parser.add_argument('-v', '--verbose', action='store_true', dest='verbose')
    pre_args, _ = pre_parser.parse_known_args(argv)

    try:
        config = load_config(pre_args.config)
    except (FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))

    if config:
        logger.debug(f"Loaded config: {config}")
        parser.set_defaults(**config)

    args = parser.parse_args(argv)

    # Configure logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.debug(f"Arguments: {args}")

    try:
        # Validate arguments
        validate_arguments(args)

        # Setup matplotlib styling
        if args.dark:
            setup_dark_background()

        # Trajectory mode: compute field along path, no plot
        if args.trajectory:
            run_trajectory(args)
            return 0

        # Plot based on planet selection
        if args.planet.lower() == 'all':
            plot_all_planets(args)
        else:
            plot_single_planet(args)

        # Handle output
        if args.save_path:
            save_plot(args.save_path)
        else:
            plt.show()

        return 0

    except ValueError as e:
        logger.error(f"Invalid argument: {e}")
        return 1
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        if args.verbose:
            logger.exception("Full traceback:")
        return 1
    finally:
        plt.close('all')


if __name__ == '__main__':
    sys.exit(main())
