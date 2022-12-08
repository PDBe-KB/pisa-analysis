import logging
import os
import os.path
#from lxml import etree 
import xml.etree.ElementTree as ET


from gemmi import cif
from pandas import DataFrame

def parse_xml_file(xml_file):
    """
    Simple helper function to read an XML file
    :param xml_file: type str - Path to XML file
    :return: parsed XML object
    """
    root = None

    xml_exists = os.path.exists(xml_file)

    if xml_exists:

        logging.debug("parsing: {}".format(xml_file))
        tree = ET.parse(xml_file)
        root = tree.getroot()

    return root


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
        int_lab_seqnum, int_resname, pdb_id, updated_cif_block
):
    """
    Function reads atom's UniProt accession numbers and sequence numbers
    from updated cif file

    :param int_lab_seqnum: type str - atom label sequence id
    :param int_resname: type str - atom residue name
    :param pdb_id: type str - PDB id
    :param updated_cif_block: parsed updated cif data
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
            
            if (
                    labseqid == int_lab_seqnum.strip()
                    and resname == int_resname.strip()
            ):
            
                n += 1
                #if labseqid != ".":
                if dbacc is not None and dbnum is not None:
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
                break
            
        
        
        if n == 0:  # If residue was not found in updated cif file, return message
            """
            label_seq_ids = updated_cif_block.find_values("_atom_site.label_seq_id")
            res_names = updated_cif_block.find_values("_atom_site.label_comp_id")
            db_acc = updated_cif_block.find_values("_atom_site.pdbx_sifts_xref_db_acc")
            db_num = updated_cif_block.find_values("_atom_site.pdbx_sifts_xref_db_num")
            
            for (labseqid, resname, dbacc, dbnum) in zip(
                    label_seq_ids,res_names, db_acc, db_num
            ):

                if (
                        labseqid == int_lab_seqnum.strip()
                        and resname == int_resname.strip()                                                                                                
                ):                                                                                                                                           
                    if dbacc is not None and dbnum is not None:

                        unp_acc = dbacc
                        unp_num = dbnum
                        #print("GDL:there are residues not mapped in fast category")
                        
                    break
            """
            if dbacc is None and dbnum is None:
                logging.debug("residue not found in updated cif file:")                                                             
                logging.debug(

                    "name {},label_seq_id {},seq_num {}, residue {}".format(                                                                             int_atname, int_lab_seqnum, int_seqnum, int_resname                                                                          )                                                                                                         
                )
                
            return unp_acc, unp_num


    
def read_cov_info(r1,r2,df: DataFrame):

    Probability=None
    Score=None
    
    if int(r1) <= int(r2) :
        unp_res1=r1
        unp_res2=r2
    if int(r2) < int(r1) :
        unp_res1=r2
        unp_res2=r1
        
    if df is not None:
        r = df[
            (df["Residue A"] == unp_res1)
            & (df["Residue B"] == unp_res2)
        ]

        Probability=list(df['Probability'])[0]
        Score=list(df['Score'])[0]

    return Probability,Score
