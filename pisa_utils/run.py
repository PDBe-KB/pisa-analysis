import argparse
import logging
import os

from pisa_utils.analyze import AnalysePisa
from pisa_utils.run_pisa import run_pisalite


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input_cif", help="input CIF file")
    parser.add_argument("--pdb_id", help="PDB ID", type=str, required=True)
    parser.add_argument("--assembly_id", help="Assembly ID", type=str, required=True)
    parser.add_argument(
        "-o", "--output_json", help="output directory for JSON and XMLs", required=True
    )
    parser.add_argument("--input_updated_cif", help="updated cif file")
    parser.add_argument("--pisa_binary", help="pisa binary path", type=str)
    parser.add_argument("--output_xml", help="output for xml files", type=str)
    parser.add_argument(
        "--force", help="always recalculate with pisa-lite", action="store_true"
    )

    args = parser.parse_args()

    logging.getLogger()

    input_cif_file = args.input_cif

    interfaces_xml_file = os.path.join(args.output_xml, "interfaces.xml")
    assembly_xml_file = os.path.join(args.output_xml, "assembly.xml")

    if (
        args.force
        or not os.path.exists(interfaces_xml_file)
        or not os.path.exists(assembly_xml_file)
    ):

        run_pisalite(input_cif=input_cif_file, output_xml=args.output_xml)

    os.makedirs(args.output_xml, exist_ok=True)

    ap = AnalysePisa(
        pdb_id=args.pdb_id,
        assembly_id=args.assembly_id,
        output_json=args.output_json,
        output_xml=args.output_xml,
        input_cif=args.input_cif,
        input_updated_cif=args.input_updated_cif,
    )

    ap.create_assem_interfaces_dict()
    ap.create_assembly_dict()


if "__main__" in __name__:
    main()
