import argparse
import logging

from pisa_utils.analyze import AnalysePisa
from pisa_utils.run_pisa import run_pisalite


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--config", help="pisa config file", required=True)
    parser.add_argument(
        "-i", "--input_dir", help="input directory with CIF directory/file"
    )
    parser.add_argument("--pdb_id", help="PDB ID", type=str, required=True)
    parser.add_argument("--assembly_id", help="Assembly ID", type=str, required=True)
    parser.add_argument(
        "-o", "--output-dir", help="output directory for JSON and XMLs", required=True
    )
    parser.add_argument("--input_updated_cif", help="updated cif file")
    parser.add_argument("--pisa_binary", help="pisa binary path", type=str)
    parser.add_argument("--result_json", help="output json file name", type=str)
    parser.add_argument(
        "--input_cif_file",
        help="assembly cif file with a seleceted name different from default",
        type=str,
    )
    parser.add_argument(
        "--force", help="always recalculate with pisa-lite", action="store_true"
    )

    args = parser.parse_args()

    logging.getLogger()

    run_pisalite(
        session_name=args.pdb_id,
        input_cif=args.input_cif_file,
        cfg_input=args.config,
        output_dir=args.output_dir,
        pisa_binary=args.pisa_binary,
    )

    # TODO: Remove obsolete arguments
    ap = AnalysePisa(
        pdb_id=args.pdb_id,
        assembly_id=args.assembly_id,
        output_dir=args.output_dir,
        result_json_file=args.result_json,
        input_dir=args.input_dir,
        input_updated_cif=args.input_updated_cif,
        input_cif_file=args.input_cif_file,
    )
    ap.process_pisa_xml()
    ap.set_results()
    ap.save_to_json()


if "__main__" in __name__:
    main()
