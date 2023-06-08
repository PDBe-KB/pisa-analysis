import os
from unittest import TestCase
from unittest.mock import patch

from pisa_utils.covariations_int import get_cov_interfaces, save_covariation_data, read_cov_info, find_element

import csv
import pandas as pd
from pandas import DataFrame
import tempfile

class TestCovariation(TestCase):


    def test_read_cov_info(self):
        """                                                                                                            
        Test that function return covariation info for two residues
                       
        """
        r1 = 24
        r2 = 52

        expected_result=(0.516928,-0.0045682)

        input_cov = os.path.join(".", "tests", "data", "F5HCH8_cov.csv")
        
        if input_cov is not None:
            covs = pd.read_csv(input_cov)
        else:
            covs = None

        a,b = read_cov_info(r1,r2,covs)

        print("GDL testing", a, b)
        self.assertEqual(read_cov_info(r1,r2,covs), expected_result)
   

    def test_read_cov_info_unpnum_none(self):
        """                                                                                                                  
        Test that function returns None when residues are not found                                                           
                                                                                                                             
        """
        r1 = 88
        r2 = 59

        expected_result=(None,None)

        input_cov = os.path.join(".", "tests", "data", "F5HCH8_cov.csv")

        if input_cov is not None:
            covs = pd.read_csv(input_cov)
        else:
            covs = None

        self.assertEqual(read_cov_info(r1,r2,covs), expected_result)


    def test_read_cov_info_df_none(self):
        """                                                                                                                                          
        Test that function returns none when data frame is empty                                                                                   
                                                                                                                                                     
        """
        r1 = 88
        r2 = 59

        expected_result=(None,None)

        input_cov = None

        if input_cov is not None:
            covs = pd.read_csv(input_cov)
        else:
            covs = None

        self.assertEqual(read_cov_info(r1,r2,covs), expected_result)

    @patch("pisa_utils.covariations_int.read_cov_info")
    def test_get_cov_interfaces(self,mock):

        """
        Test if function to overlay covariation data on interfaces returns correct list 

        """

        mock.side_effect=[(0.539227,-0.0045671), (None,None), (None,None), (None,None)]
        
        input_cov = [os.path.join(".","tests", "data", "F5HCH8_cov.csv")]
        input_json = os.path.join(".", "tests", "data", "input_json.json")
        accs_nums_list = ['F5HCH8']

        expected_result= ([['59', 'F5HCH8', 'CYS', '54', 'F5HCH8', 'CYS', 'hydrogen_bonds', '1', -0.0045671, 0.539227]],"7t4q","1")
        
        self.assertEqual(get_cov_interfaces(input_json, input_cov,accs_nums_list), expected_result)

    def test_save_covariation_data(self):
        """
        Test if function creates dataframe  with covariation data
        
        """
        
        covariation_pairs=[['59', 'F5HCH8', 'CYS', '54', 'F5HCH8', 'CYS', 'hydrogen_bonds', '1', -0.0045671, 0.539227]]
        pdb_id = "7t4q"
        assembly_id ='1'
        with tempfile.TemporaryDirectory() as output:
            expected_result = pd.DataFrame(
                covariation_pairs, columns=
                [
                    "uniprot_residue_index_a","uniprot_accession_a","residue_label_a",
                    "uniprot_residue_index_b","uniprot_accession_b","residue_label_b",
                    "contact","interface", "covariation_score", "covariation_probability"
                ]
            )
        

            result=save_covariation_data(covariation_pairs,pdb_id,assembly_id,output)
            pd.testing.assert_frame_equal(result,expected_result)

    def test_save_covariation_data_invalid_data(self):

        """
        Test function raises an error if input data is invalid
        """

        covariation_pairs=[['59', 'F5HCH8', 'CYS', '54', 'F5HCH8', 'CYS', 'hydrogen_bonds', -0.0045671, 0.539227]]

        pdb_id = "7t4q"
        assembly_id='1'

        with tempfile.TemporaryDirectory() as output:

            #save_covariation_data(covariation_pairs,pdb_id,output)
            
            self.assertRaises(TypeError, save_covariation_data,covariation_pairs)
    

    def test_save_covariation_data_return_none(self):

        """
        Test if function returns None if there is no data to save
        
        """

        covariation_pairs= []

        pdb_id = "7t4q"
        assembly_id="1"

        expected_result = None
        
        with tempfile.TemporaryDirectory() as output:
            
            result=save_covariation_data(covariation_pairs,pdb_id,assembly_id,output)

            self.assertEqual(result,expected_result)

    def test_find_element(self):

        my_array = [10, 20, 30, 40, 50]
        my_variable = 30
        expected_result = 2

        result = find_element(my_array,my_variable)
        self.assertEqual(result,expected_result)
        
