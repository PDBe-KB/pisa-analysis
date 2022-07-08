import argparse
import logging
from pisa_utils.analyse import AnalysePisa

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--config", help="pisa config file",
                        required=True)
    parser.add_argument("-i", "--input_dir",
                        help="input directory with CIF directory/file")
    parser.add_argument("--pdb_id", help="PDB ID", type=str, required=True)
    parser.add_argument("--assembly_id", help="Assembly ID", type=str,
                        required=True)
    parser.add_argument("-o", "--output-dir",
                        help="output directory for JSON and XMLs",
                        required=True)
    parser.add_argument("--input_updated_cif", help="updated cif file")
    parser.add_argument('--pisa_binary', help='pisa binary path', type=str)
    parser.add_argument('--result_json', help='output json file name', type=str)
    parser.add_argument('--input_cif_file',
                        help='assembly cif file with a seleceted name different from default',
                        type=str)
    parser.add_argument('--force', help='always recalculate with pisa-lite',
                        action='store_true')

    args = parser.parse_args()

    logging.getLogger()

    ap = AnalysePisa(pdbid_id=args.pdb_id, assembly_id=args.assembly_id,
                     pisa_config=args.config, output_dir=args.output_dir,
                     force=args.force, result_json_file=args.result_json,
                     input_dir=args.input_dir,
                     input_updated_cif=args.input_updated_cif,
                     input_cif_file=args.input_cif_file)
    ap.run_process()


if '__main__' in __name__:
    main()