import os
from lxml import etree
#import xml.etree.ElementTree as ET
from unittest import TestCase
from unittest.mock import patch

from gemmi import cif

from pisa_utils.analyze import AnalysePisa


class TestAnalyzePisa(TestCase):
    @patch("pisa_utils.analyze.get_bond_dict", side_effect=[1, 2, 3, 4, 5])
    @patch("pisa_utils.analyze.get_molecules_dict")
    #@patch("pisa_utils.analyze.get_assembly_dict")
    @patch(
        "pisa_utils.analyze.parse_xml_file",
        side_effect=[
            etree.parse(os.path.join(".", "tests", "data", "assembly.xml")).getroot(),
            etree.parse(os.path.join(".", "tests", "data", "interfaces.xml")).getroot(),
        ],
    )
    def test_create_interfaces_dict(self, bonds, molecules, xmls):
        """
        Test if the correct JSON is generated
        :param bonds: mocked bonds data
        :param molecules: mocked molecules data
        :param xmls: mocked XML data
        """
        molecules.return_value = [
            [
                {
                    "molecule_id": "1",
                    "molecule_class": "Protein",
                    "chain_id": "A-2",
                    "residue_label_comp_ids": ["ALA"],
                    "residue_seq_ids": ["2"],
                    "residue_label_seq_ids": ["5"],
                    "residue_ins_codes": [None],
                    "residue_bonds": [None],
                    "solvation_energies": [0],
                    "accessible_surface_areas": [158.39],
                    "buried_surface_areas": [0],
                }
            ],
            2,
            False,
        ]
        """
        assem_result.return_value = {
                 'assembly_id': '1',
                 'assembly_size': '6',
                 'assembly_score': '',
                 'assembly_mmsize': '2',
                 'assembly_diss_energy': 15.61,
                 'assembly_asa': 19395.3,
                 'assembly_bsa': 3514.17,
                 'assembly_entropy': 12.98,
                 'assembly_diss_area': 1427.5,
                 'assembly_int_energy': -35.28,
                 'assembly_formula': 'A(2)a(2)b(2)',
                 'assembly_composition': 'A-2A[NA](2)[GOL](2)',
                 'assem_id': '1',
                 'assembly_n_uc': '0',
                 'assembly_n_diss': '2',
                 'assembly_sym_num': '2',
                 'assembly_R350': '',
        }
        """
        
        path_to_updated_cif = os.path.join(".", "tests", "data","6nxr_updated.cif")

        ap = AnalysePisa(
            pdb_id="6nxr",
            assembly_id="1",
            output_json=os.path.join("tests", "data"),
            xmls_dir=os.path.join("tests", "data"),
            input_updated_cif=path_to_updated_cif
        )


        expected = {
            'PISA':{
                'pdb_id': '6nxr',
                'assembly_id': '1',
                'pisa_version': '2.0',
                'assembly': {
                    'mmsize': '2',
                    'dissociation_energy': 15.61,
                    'accessible_surface_area': 19395.3,
                    'buried_surface_area': 3514.17,
                    'entropy': 12.98,
                    'dissociation_area': 1427.5,
                    'solvation_energy_gain': -35.28,
                    'formula': 'A(2)a(2)b(2)',
                    'composition': 'A-2A[NA](2)[GOL](2)',
                    'interface_count': 1,
                    'interfaces':[
                        {
                            'interface_id':'1',
                            'interface_area': 1427.5,
                            'solvation_energy': -18.22,
                            'stabilization_energy': -28.59,
                            'p_value': 0.095,
                            'number_interface_residues': 2,
                            'number_hydrogen_bonds': 20,
                            'number_covalent_bonds': 0,
                            'number_disulfide_bonds': 0,
                            'number_salt_bridges': 4,
                            'number_other_bonds': 0,
                            'hydrogen_bonds': 1,
                            'salt_bridges': 2,
                            'disulfide_bonds': 4,
                            'covalent_bonds': 3,
                            'other_bonds': 5,
                            'molecules': [
                                {
                                    'molecule_id': '1',
                                    'molecule_class': 'Protein',
                                    'chain_id': 'A-2',
                                    'residue_label_comp_ids': ['ALA'],
                                    'residue_seq_ids': ['2'],
                                    'residue_label_seq_ids': ['5'],
                                    'residue_ins_codes': [None],
                                    'residue_bonds': [None],
                                    'solvation_energies': [0],
                                    'accessible_surface_areas': [158.39],
                                    'buried_surface_areas': [0]
                                }
                            ]
                        }
                    ]
                }
            }
        }

        ap.create_assem_interfaces_dict()
        self.assertEqual(ap.results, expected)


            
    
    def test_create_assembly_dict(self):

        path_to_updated_cif = os.path.join(".", "tests", "data","6nxr_updated.cif")
        ap = AnalysePisa(
            pdb_id="6nxr",
            assembly_id="1",
            output_json=os.path.join("tests", "data"),
            xmls_dir=os.path.join("tests", "data"),
            input_updated_cif=path_to_updated_cif
        )
    
        result=ap.create_assembly_dict()
    

        expected = {
            'PISA': {
                'pdb_id': '6nxr',
                'assembly_id': '1',
                'pisa_version': '2.0',
                'assembly': {
                    'id': '1',
                    'size': '6',
                    'score': '',
                    'macromolecular_size': '2',
                    'dissociation_energy': 15.61,
                    'accessible_surface_area': 19395.3,
                    'buried_surface_area': 3514.17,
                    'entropy': 12.98,
                    'dissociation_area': 1427.5,
                    'solvation_energy_gain': -35.28,
                    'number_of_uc': '0',
                    'number_of_dissociated_elements': '2',
                    'symmetry_number': '2',
                    'formula': 'A(2)a(2)b(2)',
                    'composition': 'A-2A[NA](2)[GOL](2)',
                    'R350': ''
                }
            }
        }
        
        
        self.assertEqual(result, expected)
    
    
