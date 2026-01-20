import logging
import os

from pisa_utils.analyze import AnalysePisa
from pisa_utils.constants import SUBDIR_EXTENDED_DATA
from pisa_utils.models.models import LigandProcessingMode
from pisa_utils.parsers import (
    CompileInterfaceSummaryJSON,
    ConvertAssemblyListToJSON,
    ConvertAssemblyXMLToJSON,
    ConvertComponentsListToJSON,
    ConvertInterfaceListToJSON,
    ConvertInterfaceXMLToJSONs,
)
from pisa_utils.run_pisa import run_pisa_service, run_pisalite
from pisa_utils.utils import collect_base_args, validate_args


def main():
    """CLI entry point for running PISA analysis on a given assembly CIF file"""

    parser = collect_base_args()
    parser.add_argument("--pdb_id", help="PDB ID", required=True)
    parser.add_argument("--assembly_id", help="Assembly ID", required=True)
    parser.add_argument("--input_updated_cif", help="updated cif path/file")

    args = parser.parse_args()
    args = validate_args(args)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)

    interfaces_xml_file = os.path.join(args.output_xml, "interfaces.xml")
    assembly_xml_file = os.path.join(args.output_xml, "assembly.xml")
    assembly_json = os.path.join(
        args.output_json, f"{args.pdb_id}-assembly{args.assembly_id}.json"
    )
    interface_json = os.path.join(
        args.output_json, f"{args.pdb_id}-assembly{args.assembly_id}-interfaces.json"
    )

    logging.info("input cif: {}".format(args.input_cif))
    logging.info("output xml: {}".format(args.output_xml))
    logging.info("output json: {}".format(args.output_json))

    if (
        args.force
        or not os.path.exists(interfaces_xml_file)
        or not os.path.exists(assembly_xml_file)
    ):

        run_pisalite(
            input_cif=args.input_cif,
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

    interfaces = ap.interfaces_xml_to_json(
        assembly_xml_file, interfaces_xml_file, interface_json
    )
    ap.assembly_xml_to_json(assembly_xml_file, assembly_json, interfaces)


def service():
    """
    CLI entry point for running PISA with additional options. Intended for internal
    use at the PDBe.
    """

    parser = collect_base_args()
    parser.add_argument(
        "--as_is", help="run pisa in asis mode", action="store_true", default=False
    )
    parser.add_argument(
        "--lig_proc_mode",
        help="ligand processing mode for pisa",
        default=LigandProcessingMode.AUTO.value,
        choices=[mode.value for mode in LigandProcessingMode],
    )
    parser.add_argument(
        "--exclude_ligands",
        help="list of ligand CCD (three-letter codes) to exclude",
        nargs="+",
        default=None,
    )
    parser.add_argument(
        "--run_all_pisa_commands",
        help=(
            "parse full results from pisa XMLs, including all assemblies, interfaces, "
            "and monomers"
        ),
        action="store_true",
        default=True,
    )
    args = parser.parse_args()
    args = validate_args(args)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)

    logging.info("input cif: {}".format(args.input_cif))
    logging.info("output xml: {}".format(args.output_xml))
    logging.info("output json: {}".format(args.output_json))

    os.makedirs(args.output_xml, exist_ok=True)

    (
        assembly_xml_file,
        interfaces_xml_file,
        assemblies_extended_file,
        interfaces_extended_file,
        monomers_extended_file,
    ) = run_pisa_service(
        input_cif=args.input_cif,
        xml_output_dir=args.output_xml,
        pisa_binary=args.pisa_binary,
        pisa_setup_dir=args.pisa_setup_dir,
        asis=args.as_is,
        lig_proc_mode=args.lig_proc_mode,
        exclude_ligands=args.exclude_ligands,
        extended_data=args.run_all_pisa_commands,
    )

    if args.run_all_pisa_commands:

        data_parser = (
            ConvertInterfaceXMLToJSONs(
                path_xml=interfaces_xml_file,
                path_jsons=os.path.join(args.output_json, "interfaces"),
                path_structure_file=args.input_cif,
            ),
            ConvertAssemblyXMLToJSON(
                path_xml=assembly_xml_file,
                path_json=os.path.join(args.output_json, "assemblies.json"),
                path_structure_file=args.input_cif,
                path_interface_jsons=os.path.join(args.output_json, "interfaces"),
            ),
            CompileInterfaceSummaryJSON(
                path_interface_jsons=os.path.join(args.output_json, "interfaces"),
                path_assembly_json=os.path.join(args.output_json, "assemblies.json"),
                path_output_json=os.path.join(
                    args.output_json, "interface_summary.json"
                ),
            ),
            ConvertAssemblyListToJSON(
                path_txt=assemblies_extended_file,
                path_json=os.path.join(args.output_json, SUBDIR_EXTENDED_DATA),
            ),
            ConvertInterfaceListToJSON(
                path_txt=interfaces_extended_file,
                path_json=os.path.join(args.output_json, SUBDIR_EXTENDED_DATA),
            ),
            ConvertComponentsListToJSON(
                path_txt=monomers_extended_file,
                path_json=os.path.join(args.output_json, SUBDIR_EXTENDED_DATA),
            ),
        )

        for parser in data_parser:
            logging.info(f"Parsing and converting: {parser}")
            parser.parse()


if __name__ == "__main__":
    main()
