import json
import logging
import os
import tempfile
from unittest import TestCase
from unittest.mock import patch

from pisa_utils.analyze import AnalysePisa

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.DEBUG)


mol_dict = (
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
            "solvation_energies": [0.0],
            "accessible_surface_areas": [158.39],
            "buried_surface_areas": [0],
        }
    ],
    2,
    False,
)


class TestAnalyzePisa(TestCase):
    @patch("pisa_utils.analyze.get_bond_dict", side_effect=[1, 2, 3, 4, 5])
    @patch("pisa_utils.analyze.get_molecules_dict", return_value=mol_dict)
    def test_create_interfaces_dict(self, bonds, molecules):
        """
        Test if the correct JSON is generated
        :param bonds: mocked bonds data
        :param molecules: mocked molecules data
        :param xmls: mocked XML data
        """

        path_to_updated_cif = os.path.join(".", "tests", "data", "6nxr_updated.cif")

        ap = AnalysePisa(
            pdb_id="6nxr",
            assembly_id="1",
            input_updated_cif=path_to_updated_cif,
        )

        xmls_dir = os.path.join("tests", "data")

        assembly_xml = os.path.join(xmls_dir, "non-existent.xml")
        interface_xml = os.path.join(xmls_dir, "non-existent-interface.xml")

        with tempfile.NamedTemporaryFile(mode="wt") as f:
            with self.assertRaises(FileNotFoundError):
                ap.interfaces_xml_to_json(assembly_xml, interface_xml, f.name)

        assembly_xml = os.path.join(xmls_dir, "assembly.xml")
        interface_xml = os.path.join(xmls_dir, "non-existent.xml")

        with tempfile.NamedTemporaryFile(mode="wt") as f:
            with self.assertRaises(FileNotFoundError):
                ap.interfaces_xml_to_json(assembly_xml, interface_xml, f.name)

        assembly_xml = os.path.join(xmls_dir, "non-existent.xml")
        interface_xml = os.path.join(xmls_dir, "interfaces.xml")
        with tempfile.NamedTemporaryFile(mode="wt") as f:
            with self.assertRaises(FileNotFoundError):
                ap.interfaces_xml_to_json(assembly_xml, interface_xml, f.name)

        assembly_xml = os.path.join(xmls_dir, "assembly.xml")
        interface_xml = os.path.join(xmls_dir, "interfaces.xml")
        with tempfile.NamedTemporaryFile(mode="wt") as f:
            expected = json.load(
                open(os.path.join("tests", "data", "6nxr-assembly1-interfaces.json"))
            )
            expected = {
                "6nxr": {
                    "assembly_id": "1",
                    "pisa_version": "2.0",
                    "assembly": {
                        "mmsize": "2",
                        "dissociation_energy": 15.61,
                        "accessible_surface_area": 19395.3,
                        "buried_surface_area": 3514.17,
                        "entropy": 12.98,
                        "dissociation_area": 1427.5,
                        "solvation_energy_gain": -35.28,
                        "formula": "A(2)a(2)b(2)",
                        "composition": "A-2A[NA](2)[GOL](2)",
                        "interface_count": 1,
                        "interfaces": [
                            {
                                "interface_id": "1",
                                "interface_area": 1427.5,
                                "solvation_energy": -18.22,
                                "stabilization_energy": -28.59,
                                "p_value": 0.095,
                                "number_interface_residues": 2,
                                "number_hydrogen_bonds": 20,
                                "number_covalent_bonds": 0,
                                "number_disulfide_bonds": 0,
                                "number_salt_bridges": 4,
                                "number_other_bonds": 0,
                                "hydrogen_bonds": 1,
                                "salt_bridges": 2,
                                "disulfide_bonds": 4,
                                "covalent_bonds": 3,
                                "other_bonds": 5,
                                "molecules": [
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
                            }
                        ],
                    },
                }
            }
            ap.interfaces_xml_to_json(assembly_xml, interface_xml, f.name)
            molecules.assert_called()
            bonds.assert_called()
            result = json.load(open(f.name))
            self.maxDiff = None
            self.assertEqual(result, expected)

    def test_create_assembly_dict(self):

        path_to_updated_cif = os.path.join(".", "tests", "data", "6nxr_updated.cif")
        ap = AnalysePisa(
            pdb_id="6nxr", assembly_id="1", input_updated_cif=path_to_updated_cif
        )

        expected = {
            "6nxr": {
                "assembly_id": "1",
                "pisa_version": "2.0",
                "assembly": {
                    "id": "1",
                    "size": "6",
                    "interface_count": 1,
                    "score": "",
                    "macromolecular_size": "2",
                    "dissociation_energy": 15.61,
                    "accessible_surface_area": 19395.3,
                    "buried_surface_area": 3514.17,
                    "entropy": 12.98,
                    "dissociation_area": 1427.5,
                    "solvation_energy_gain": -35.28,
                    "number_of_uc": "0",
                    "number_of_dissociated_elements": "2",
                    "symmetry_number": "2",
                    "formula": "A(2)a(2)b(2)",
                    "composition": "A-2A[NA](2)[GOL](2)",
                    "R350": "",
                },
            }
        }
        with tempfile.NamedTemporaryFile(mode="wt") as f:
            interfaces = 1
            assembly_xml = os.path.join("tests", "data", "assembly.xml")
            ap.assembly_xml_to_json(assembly_xml, f.name, interfaces)
            self.assertEqual(json.load(open(f.name, "r")), expected)
