import argparse
import json
import logging
import os
from time import time
from gemmi import cif
import csv

import pandas as pd
from pandas import DataFrame

#from pisa_utils.utils import read_cov_info


def get_cov_interfaces(input_json, input_cov):

    """
    Add covariation data to interface contacts

    Creates csv file with covariation info for interfaces' contacts 
                                                                                                   
    Args:                                                                                                        
        input_json (str): json file with assembly interfaces data (obtained with pisa-analysis)

        input_cov (str) : csv file with covariation pairs (probability >= 0.5) for uniprot acc

        output (str) : Path to output file containing interfaces contacts with covariation signal 

    """                                                                                              
    # read interfaces info
    
    fp = open(input_json, "r")
    data = json.load(fp)
    fp.close()
    pisa=list(data.keys())[0]
    data = data[pisa]
    pdb_id = data['pdb_id']
    assembly_id = data['assembly_id']
    assembly_data = data['assembly']
    
    interfaces = assembly_data['interfaces']
    
    #read covariation pairs for uniprot accession                                                            
    if input_cov is not None:
        covs = pd.read_csv(input_cov)
    else:
        covs = None

    covariation_pairs=[]
       
    for interface in interfaces:
        interface_id=interface.get('interface_id')
        
        for key, prop in interface.items():
            if key in [
                    "hydrogen_bonds",
                    "salt_bridges",
                    "disulfide_bonds",
                    "covalent_bonds",
                    "other_bonds"
            ]:
                unp_nums_1=prop["atom_site_1_unp_nums"]
                unp_nums_2=prop["atom_site_2_unp_nums"]
                unp_accs_1=prop["atom_site_2_unp_accs"]
                unp_accs_2=prop["atom_site_2_unp_accs"]
                residues_1=prop["atom_site_1_residues"]
                residues_2=prop["atom_site_2_residues"]
                
                for ( unp_num_1, unp_num_2, unp_acc_1,
                      unp_acc_2,residue_1, residue_2 ) in zip(
                      unp_nums_1,unp_nums_2,unp_accs_1,
                      unp_accs_2,residues_1,residues_2):

                     if (unp_num_1 is not None) and (unp_num_2 is not None):

                        cov_probability, cov_score = read_cov_info(unp_num_1,unp_num_2,covs)
                        
                        if (cov_probability is not None) and (cov_score is not None):
                            covariation_pairs.append(
                                [
                                    unp_num_1,unp_acc_1,residue_1,
                                    unp_num_2,unp_acc_2,residue_2,
                                    key,interface_id,cov_score,
                                    cov_probability
                                ]
                            )
                        
     
    return covariation_pairs, pdb_id

def save_covariation_data(covariation_pairs,pdb_id,output_dir):
    """
    Creates output csv file with interfaces contacts with covariation signal

    Args:
    
    covariation_pairs (list) : List of assembly interfaces with covariation signal

    pdb_id : Entry pdb id

    output_dir (str) : Path to output directory to save file containing interfaces contacts with covariation signal
    
    """
    if len(covariation_pairs) != 0 :

        try:
        
            df = pd.DataFrame(
                covariation_pairs, columns=
                [
                    "uniprot_residue_index_a","uniprot_accession_a","residue_label_a",
                    "uniprot_residue_index_b","uniprot_accession_b","residue_label_b",
                    "contact","interface", "covariation_score", "covariation_probability"
                ]
            )
            
            output_csv = os.path.join(output_dir,"{}_interfaces_cov.csv".format(pdb_id))
            df.to_csv(output_csv, index=False)

            return df 


        except Exception as e:
            logging.error(
                "Invalid data frame for covariation pairs: probably wrong fields in the data "
            )
            logging.error(e)
    else:
        logging.info(f"No covariation data found to save ' contacts")

        return None 


def read_cov_info(r1,r2,df: DataFrame):

    """
    Outputs covariation data for residues r1 and r2
    
    Args:

    r1 : uniprot sequence number for residue 1

    r2 : uniprot sequence number for residue 2

    df: DataFrame : data containing covariation information for residue pairs 
    
    """
    
    Probability=None
    Score=None

    if int(r1) <= int(r2) :
        unp_res1=int(r1)
        unp_res2=int(r2)
    if int(r2) < int(r1) :
        unp_res1=int(r2)
        unp_res2=int(r1)

    if df is not None:

        n=0
        index=None
        for i,j in zip(df['uniprot_residue_index_a'],df['uniprot_residue_index_b']):
            if i==unp_res1 and j==unp_res2:
                index=n
            n=n+1
        if index is not None:
            
            Probability = df['covariation_probability'][index]
            Score = df['covariation_score'][index]

    
    return Probability,Score


def main():
    """Application entry point"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input_json",
        help="Input json file with interfaces information",
        required=True,
    )
    parser.add_argument(
        "--input_cov",
        help="Input covariation pairs file for uniprot accession",
        required=True,
    )
    parser.add_argument(
        "-o", "--output_csv", help="output directory for csv file", required=True
    )

    args = parser.parse_args()
    
    covariation_pairs, pdb_id = get_cov_interfaces(args.input_json, args.input_cov)

    save_covariation_data(covariation_pairs,pdb_id,args.output_csv)


    logging.info("We are done here.")

if __name__ == "__main__":
    main()

