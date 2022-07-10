import os
import xml.etree.ElementTree as ET
from unittest import TestCase

from pisa_utils.dictionaries import get_molecules_dict


class TestDictionaries(TestCase):
    def test_get_molecules_dict(self):
        """
        Test that the function returns a correct dictionary
        """
        molecule = ET.parse(
            os.path.join(".", "tests", "data", "mocks", "interfaces.xml")
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
            os.path.join(".", "tests", "data", "mocks", "interfaces.xml")
        ).getroot()
        molecules = [molecule]
        self.assertEqual(get_molecules_dict(molecules)[1], 1)

    def test_get_molecules_invalid(self):
        """
        Test that the function returns invalid when only 1 interface
        """
        molecule = ET.parse(
            os.path.join(".", "tests", "data", "mocks", "interfaces.xml")
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
