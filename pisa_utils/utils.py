import logging
import os
import xml.etree.ElementTree as ET

from gemmi import cif


def parse_xml_file(xml_file):
    """
    Simple helper function to read an XML file
    :param xml_file: type str - Path to XML file
    :return: parsed XML object
    """
    root = None
    try:
        logging.debug("parsing: {}".format(xml_file))
        tree = ET.parse(xml_file)
        root = tree.getroot()
    except Exception as e:
        logging.error("invalid xml file: {}".format(xml_file))
        logging.error(e)
    return root

def create_pisa_config(dataroot,setup_dir):
    """
    Function creates a configuration file to run pisa

    :param dataroot : path to output pisa data and session
    :param setup_dir : path to setup directory of pisa-lite
    :return : None

    """
    infile=open(os.path.join(setup_dir,"pisa_cfg_tmp"))
    outputname=os.path.join(dataroot,"pisa.cfg")
    outfile=open(outputname,"w")
    for line in infile:
        line=line.replace("path_dataroot",dataroot).replace("path_to_setup",setup_dir)
        outfile.write(line)
    infile.close()
    outfile.close()

def read_uniprot_info(
    int_lab_seqnum, int_seqnum, int_atname, int_resname, pdb_id, input_updated_cif
):
    """
    Function reads atom's UniProt accession numbers and sequence numbers
    from updated cif file

    :param int_lab_seqnum: type str - atom label sequence id
    :param int_seqnum: type str - atom sequence number
    :param int_atname: type str - atom name
    :param int_resname: type str - atom residue name
    :param pdb_id: type str - PDB id
    :param input_updated_cif: type str - path to updated cif file
    :return: (atom uniprot accession number, atom uniprot sequence number)
    """
    unp_acc = None
    unp_num = None

    # updated cif file path
    path = os.path.join(input_updated_cif, "{}_updated.cif".format(pdb_id))

    # Reading UniProt acc and seq numbers in updated cif file:
    doc = cif.read_file(path)
    block = doc.sole_block()
    label_seq_id = block.find_loop("_atom_site.label_seq_id")
    atom_name = block.find_loop("_atom_site.label_atom_id")
    auth_seq_id = block.find_loop("_atom_site.auth_seq_id")
    db_name = block.find_loop("_atom_site.pdbx_sifts_xref_db_name")
    db_acc = block.find_loop("_atom_site.pdbx_sifts_xref_db_acc")
    db_num = block.find_loop("_atom_site.pdbx_sifts_xref_db_num")
    res_name = block.find_loop("_atom_site.label_comp_id")

    # counts if atom is not found in updated cif file
    n = 0

    # if sequence identifier read in pisa-lite is not available,
    # replace with a dot
    if int(int_lab_seqnum) < 0:
        int_lab_seqnum = "."

    # Search atom in updated cif file and read uniprot numbers
    for (labseqnum, seqnum, name, resname, dbname, dbacc, dbnum) in zip(
        label_seq_id, auth_seq_id, atom_name, res_name, db_name, db_acc, db_num
    ):
        if (
            labseqnum == int_lab_seqnum.strip()
            and seqnum == int_seqnum.strip()
            and name == int_atname.strip()
            and resname == int_resname.strip()
        ):
            n += 1
            if dbname == "UNP":
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
    if n == 0:  # If atom was not found in updated cif file, return message
        logging.debug("atom not found in updated cif file:")
        logging.debug(
            "name {},label_seq_id {},seq_num {}, residue {}".format(
                int_atname, int_lab_seqnum, int_seqnum, int_resname
            )
        )
        return unp_acc, unp_num
