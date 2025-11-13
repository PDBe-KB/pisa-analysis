import argparse
import logging
import os

from pisa_utils.analyze import AnalysePisa
from pisa_utils.run_pisa import run_pisalite


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input_cif",
        help="Input CIF file",
        required=True,
    )
    parser.add_argument("--pdb_id", help="PDB ID", required=True)
    parser.add_argument("--assembly_id", help="Assembly ID", required=True)
    parser.add_argument(
        "-o", "--output_json", help="output directory for JSON and XMLs", required=True
    )
    parser.add_argument("--input_updated_cif", help="updated cif path/file")
    parser.add_argument(
        "--output_xml", help="output directory for xml files", required=True
    )
    parser.add_argument("--pisa_binary", help="path to pisa binary", default="pisa")
    parser.add_argument(
        "--pisa_setup_dir",
        help="path to pisa-lite setup directory",
    )
    parser.add_argument(
        "--force", help="always recalculate with pisa-lite", action="store_true"
    )
    parser.add_argument(
        "-v", "--verbose", help="Add debug information to log", action="store_true"
    )

    args = parser.parse_args()
    
    if not args.pisa_setup_dir:
        args.pisa_setup_dir = os.environ.get("PISA_SETUP_DIR")

    if not args.pisa_setup_dir:
        raise Exception(
            "PISA_SETUP_DIR not set in environment and --pisa_setup_dir not specified"
        )
    #if not (args.pisa_setup_dir or "PISA_SETUP_DIR" in os.environ):
    #    raise Exception(
    #        "PISA_SETUP_DIR not set in environment and --pisa_setup_dir not specified"
    #    )

    if not "pisa" and not os.path.isfile(args.pisa_binary):
        raise Exception(f"pisa binary not found or is not a file: {args.pisa_binary}")

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)

    input_cif_file = args.input_cif

    interfaces_xml_file = os.path.join(args.output_xml, "interfaces.xml")
    assembly_xml_file = os.path.join(args.output_xml, "assembly.xml")
    assembly_json = os.path.join(
        args.output_json, f"{args.pdb_id}-assembly{args.assembly_id}.json"
    )
    interface_json = os.path.join(
        args.output_json, f"{args.pdb_id}-assembly{args.assembly_id}-interfaces.json"
    )

    logging.info("input cif: {}".format(input_cif_file))
    logging.info("output xml: {}".format(args.output_xml))
    logging.info("output json: {}".format(args.output_json))

    if (
        args.force
        or not os.path.exists(interfaces_xml_file)
        or not os.path.exists(assembly_xml_file)
    ):

        run_pisalite(
            input_cif=input_cif_file,
            xml_output_dir=args.output_xml,
            pisa_binary=args.pisa_binary,
            pisa_setup_dir=args.pisa_setup_dir,
        )

    os.makedirs(args.output_xml, exist_ok=True)

    ap = AnalysePisa(
        pdb_id=args.pdb_id,
        assembly_id=args.assembly_id,
        input_updated_cif=args.input_updated_cif,
    )

    
    interfaces = ap.interfaces_xml_to_json(assembly_xml_file, interfaces_xml_file, interface_json)
    ap.assembly_xml_to_json(assembly_xml_file, assembly_json,interfaces)

if "__main__" in __name__:
    main()
