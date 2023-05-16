import logging
import os
import os.path
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


def create_pisa_config(dataroot, setup_dir):
    """
    Function creates a configuration file to run pisa

    :param dataroot : path to output pisa data and session
    :param setup_dir : path to setup directory of pisa-lite
    :return : None

    """
    outputname = os.path.join(dataroot, "pisa.cfg")
    with open(os.path.join(setup_dir, "pisa_cfg_tmp")) as infile:
        with open(outputname, "w") as outfile:
            for line in infile:
                line = line.replace("path_dataroot", dataroot).replace(
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

    # if sequence identifier read in pisa-lite is not available,
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

        for (labseqid, resname, dbacc, dbnum) in zip(
            db_seq_id, res_name, db_acc, db_num
        ):

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
