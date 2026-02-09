import argparse
import logging
import os
import os.path
import shutil
import xml.etree.ElementTree as ET


def parse_xml_file(xml_file):
    """Parse XML file.

    Args:
        xml_file (str): Path to XML file.

    Returns:
        xml.etree.ElementTree.Element: Root element of XML file.

    """
    tree = ET.parse(xml_file)
    return tree.getroot()


def create_pisa_config(dataroot, setup_dir, placeholder_dataroot="path_dataroot"):
    """
    Function creates a configuration file to run pisa

    :param dataroot : path to output pisa data and session
    :param setup_dir : path to setup directory of pisa
    :return : None

    """
    outputname = os.path.join(dataroot, "pisa.cfg")
    with open(os.path.join(setup_dir, "pisa_cfg_tmp")) as infile:
        with open(outputname, "w") as outfile:
            for line in infile:
                line = line.replace(placeholder_dataroot, dataroot).replace(
                    "path_to_setup", setup_dir
                )
                outfile.write(line)
    return outputname


def read_uniprot_info(
    int_lab_seqnum, int_seqnum, int_atname, int_resname, updated_cif_block
):
    """
    Function reads atom's UniProt accession numbers and sequence numbers
    from updated cif file

    :param int_lab_seqnum: type str - atom label sequence id
    :param int_seqnum: type str - atom sequence number
    :param int_atname: type str - atom name
    :param int_resname: type str - atom residue name
    :param parsed updated cif data: type str - path to updated cif file
    :return: (atom uniprot accession number, atom uniprot sequence number)
    """

    unp_acc = None
    unp_num = None

    # counts if atom is not found in updated cif file
    n = 0

    # if sequence identifier read in pisa is not available,
    # replace with a dot
    if int(int_lab_seqnum) < 0:
        int_lab_seqnum = "."

    if not (updated_cif_block):
        return unp_acc, unp_num

    else:
        db_seq_id = updated_cif_block.find_values("_pdbx_sifts_xref_db.seq_id")
        db_acc = updated_cif_block.find_values("_pdbx_sifts_xref_db.unp_acc")
        db_num = updated_cif_block.find_values("_pdbx_sifts_xref_db.unp_num")
        res_name = updated_cif_block.find_values("_pdbx_sifts_xref_db.mon_id")

        for labseqid, resname, dbacc, dbnum in zip(db_seq_id, res_name, db_acc, db_num):
            if labseqid == int_lab_seqnum.strip() and resname == int_resname.strip():
                n += 1
                if labseqid != ".":
                    unp_acc = dbacc
                    unp_num = dbnum

                    return unp_acc, unp_num

                else:
                    logging.debug("No UNP numbers found for atom:")
                    logging.debug(
                        "name {},label_seq_id {},seq_num {}, residue {}".format(
                            int_atname, int_lab_seqnum, int_seqnum, int_resname
                        )
                    )
                    unp_acc = None
                    unp_num = None
                    return unp_acc, unp_num
        # If residue was not found in updated cif file, return message
        if not n:
            logging.debug("residue not found in updated cif file:")
            logging.debug(
                "name {},label_seq_id {},seq_num {}, residue {}".format(
                    int_atname, int_lab_seqnum, int_seqnum, int_resname
                )
            )
            return unp_acc, unp_num


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


def extract_ligand_contents(ligand_id: str) -> tuple[str, str, int]:
    """
    Extract auth_asym_id, ccd_id, and auth_seq_id from ligand identifier.

    ligand_id usually formatted as: [CCD_ID]auth_asym_id:auth_seq_id, e.g. [ANF]H:500

    :param ligand_id: Ligand identifier. e.g. [ANF]H:500
    :type ligand_id: str
    :return: auth_asym_id, ccd_id, auth_seq_id.
    :rtype: tuple[str, str, int]
    """
    chain_id_split = ligand_id.split(":")
    auth_asym_id = chain_id_split[0].split("]")[1]
    ccd_id = ligand_id.split("]")[0][1:]
    auth_seq_id = int(chain_id_split[1])

    return auth_asym_id, ccd_id, auth_seq_id


def id_is_ligand(id: str) -> bool:
    """
    Determine if an identifier corresponds to a ligand.

    Usually formatted as: [CCD_ID]auth_asym_id:auth_seq_id, e.g. [ANF]H:500

    :param id: Identifier to check. e.g. [ANF]H:500
    :type id: str
    :return: True if identifier corresponds to a ligand, False otherwise.
    :rtype: bool
    """
    return "[" in id and "]" in id and ":" in id


def is_int_or_float(value: str) -> bool:
    """
    Check if a string can be converted to an integer or a float.

    :param value: String to check.
    :type value: str
    :return: True if the string can be converted to an integer or a float, False
        otherwise.
    :rtype: bool
    """
    # Integer
    if value.isdigit():
        return True
    # Float
    try:
        float(value)
        return True
    except ValueError:
        return False
