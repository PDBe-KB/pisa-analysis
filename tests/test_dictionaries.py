import os
import xml.etree.ElementTree as ET
from unittest import TestCase
from unittest.mock import patch

from pisa_utils.dictionaries import get_bond_dict, get_molecules_dict


class TestDictionaries(TestCase):

    # Note that the mocked function is imported in dictionaries
    @patch("pisa_utils.dictionaries.read_uniprot_info")
    def test_get_bond_dict(self, mock):
        mock.return_value = ("P48491", "2")
        bond = ET.parse(
            os.path.join(".", "tests", "data", "mocks", "bonds.xml")
        ).getroot()
        bonds = [bond]
        result = get_bond_dict(bonds, "H-bond", "pdbid", "cifpath")
        expected = {
            "bond_distances": [2.99],
            "atom_site_1_chains": ["A-2"],
            "atom_site_1_residues": ["ASN"],
            "atom_site_1_label_asym_ids": ["A-2"],
            "atom_site_1_orig_label_asym_ids": ["A"],
            "atom_site_1_unp_accs": ["P48491"],
            "atom_site_1_unp_nums": ["2"],
            "atom_site_1_seq_nums": ["10"],
            "atom_site_1_label_seq_ids": ["13"],
            "atom_site_1_label_atom_ids": ["ND2"],
            "atom_site_1_inscodes": [None],
            "atom_site_2_chains": ["A"],
            "atom_site_2_residues": ["THR"],
            "atom_site_2_label_asym_ids": ["A"],
            "atom_site_2_orig_label_asym_ids": ["A"],
            "atom_site_2_unp_accs": ["P48491"],
            "atom_site_2_unp_nums": ["2"],
            "atom_site_2_seq_nums": ["76"],
            "atom_site_2_label_seq_ids": ["79"],
            "atom_site_2_label_atom_ids": ["OG1"],
            "atom_site_2_inscodes": [None],
        }

        self.assertEqual(result, expected)

    def test_get_molecules_dict(self):
        """
        Test that the function returns a correct dictionary
        """
        molecule = ET.parse(
            os.path.join(".", "tests", "data", "mocks", "molecules.xml")
        ).getroot()
        molecules = [molecule]
        expected_dictionary = [
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
        ]
        self.assertEqual(get_molecules_dict(molecules)[0][0], expected_dictionary[0])

    def test_get_molecules_counts(self):
        """
        Test that the function returns a correct count
        """
        molecule = ET.parse(
            os.path.join(".", "tests", "data", "mocks", "molecules.xml")
        ).getroot()
        molecules = [molecule]
        self.assertEqual(get_molecules_dict(molecules)[1], 1)

    def test_get_molecules_invalid(self):
        """
        Test that the function returns invalid when only 1 interface
        """
        molecule = ET.parse(
            os.path.join(".", "tests", "data", "mocks", "molecules.xml")
        ).getroot()
        molecules = [molecule]
        self.assertTrue(get_molecules_dict(molecules)[2])

    def test_get_molecules_ligand(self):
        """
        Test that the function returns invalid when ligand interface
        """
        molecule = ET.parse(
            os.path.join(".", "tests", "data", "mocks", "interfaces_ligand.xml")
        ).getroot()
        molecules = [molecule]
        self.assertTrue(get_molecules_dict(molecules)[2])
