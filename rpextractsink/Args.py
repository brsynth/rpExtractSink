from argparse import ArgumentParser

from brs_utils import add_logger_args

DEFAULTS = {"comp": "c", "cspace": "mnx4.4"}


def add_arguments(parser: ArgumentParser) -> ArgumentParser:
    # Add arguments related to the logger
    parser = add_logger_args(parser)

    parser.add_argument("input_sbml", type=str, help="input SBML file")
    parser.add_argument("output_sink", type=str, help="output sink file")
    parser.add_argument(
        "--compartment-id",
        type=str,
        default=DEFAULTS["comp"],
        help=(
            f"SBML compartment id from which to extract "
            f"the chemical species (default: {DEFAULTS['comp']})"
        ),
    )
    parser.add_argument(
        "--remove-dead-end",
        action="store_true",
        help=(
            "upon FVA evaluation, ignore chemical" "species that do not have any flux"
        ),
    )
    parser.add_argument(
        "--chemical-space",
        dest="cspace",
        type=str,
        default=DEFAULTS["cspace"],
        help=(f"Type of cache data to use (default: {DEFAULTS['cspace']})"),
    )

    parser.add_argument(
        "--cache-dir",
        dest="cache_dir",
        default=None,
        type=str,
        help="Directory where the rrCache is stored. If not provided, the default directory will be used.",
    )
    parser.add_argument(
        "--standalone", action="store_true", help=("do not get the InChI from Internet")
    )
    return parser
