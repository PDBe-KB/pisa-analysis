import argparse
import os
import os.path
import shutil


def collect_base_args() -> argparse.ArgumentParser:
    """Adds base arguments to an ArgumentParser object."""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input_cif",
        help="Input CIF file",
        required=True,
    )
    parser.add_argument(
        "-o", "--output_json", help="output directory for JSON and XMLs", required=True
    )
    parser.add_argument(
        "--output_xml", help="output directory for xml files", required=True
    )
    parser.add_argument("--pisa_binary", help="path to pisa binary", default="pisa")
    parser.add_argument(
        "--pisa_setup_dir",
        help="path to pisa setup directory",
    )
    parser.add_argument(
        "--force", help="always recalculate with pisa", action="store_true"
    )
    parser.add_argument(
        "-v", "--verbose", help="Add debug information to log", action="store_true"
    )

    return parser


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """Checks command line arguments for validity."""

    if not args.pisa_setup_dir:
        args.pisa_setup_dir = os.environ.get("PISA_SETUP_DIR")

    if not args.pisa_setup_dir:
        raise Exception(
            "PISA_SETUP_DIR not set in environment and --pisa_setup_dir not specified"
        )

    if not (os.path.isfile(args.pisa_binary) or shutil.which(args.pisa_binary)):
        raise Exception(f"pisa binary not found or is not a file: {args.pisa_binary}")

    return args
