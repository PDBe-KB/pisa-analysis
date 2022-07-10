import os
import xml.etree.ElementTree as ET
from unittest import TestCase
from unittest.mock import patch

from pisa_utils.analyze import AnalysePisa


class TestAnalyzePisa(TestCase):
    @patch("pisa_utils.analyze.get_bond_dict", side_effect=[1, 2, 3, 4, 5])
    @patch("pisa_utils.analyze.get_molecules_dict")
    @patch(
        "pisa_utils.analyze.parse_xml_file",
        side_effect=[
            ET.parse(os.path.join(".", "tests", "data", "assembly.xml")).getroot(),
            ET.parse(os.path.join(".", "tests", "data", "interfaces.xml")).getroot(),
        ],
    )
    def test_process_pisa_xml(self, bonds, molecules, xmls):
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
        ap = AnalysePisa(
            pdb_id="6nxr",
            assembly_id="1",
            output_dir=os.path.join("tests", "data"),
            result_json_file=".",
            input_dir=".",
            input_updated_cif=".",
            input_cif_file=".",
        )
        expected = {
            "assembly_status": "Ok",
            "status": "Ok",
            "num_interfaces": "5",
            "id": ["1"],
            "int_area": [1427.5],
            "interface_dicts": [
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
            "non_ligand_interface_count": 1,
            "assembly_mmsize": "2",
            "assembly_diss_energy": 15.61,
            "assembly_asa": 19395.3,
            "assembly_bsa": 3514.17,
            "assembly_entry": 12.98,
            "assembly_diss_area": 1427.5,
            "assembly_int_energy": -35.28,
            "assembly_formula": "A(2)a(2)b(2)",
            "assembly_composition": "A-2A[NA](2)[GOL](2)",
        }
        ap.process_pisa_xml()
        self.assertEqual(ap.interfaces_results, expected)

    def test_set_results(self):
        ap = AnalysePisa(
            pdb_id="6nxr",
            assembly_id="1",
            output_dir=os.path.join("tests", "data"),
            result_json_file=".",
            input_dir=".",
            input_updated_cif=".",
            input_cif_file=".",
        )
        ap.interfaces_results = {
            "assembly_status": "Ok",
            "status": "Ok",
            "num_interfaces": "5",
            "id": ["1"],
            "int_area": [1427.5],
            "interface_dicts": [
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
            "non_ligand_interface_count": 1,
            "assembly_mmsize": "2",
            "assembly_diss_energy": 15.61,
            "assembly_asa": 19395.3,
            "assembly_bsa": 3514.17,
            "assembly_entry": 12.98,
            "assembly_diss_area": 1427.5,
            "assembly_int_energy": -35.28,
            "assembly_formula": "A(2)a(2)b(2)",
            "assembly_composition": "A-2A[NA](2)[GOL](2)",
        }
        ap.set_results()
        expected = {
            "PISA": {
                "pdb_id": "6nxr",
                "assembly_id": "1",
                "pisa_version": "2.0",
                "assembly": {
                    "mmsize": "2",
                    "dissociation_energy": 15.61,
                    "accessible_surface_area": 19395.3,
                    "buried_surface_area": 3514.17,
                    "entropy": None,
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
        self.assertEqual(ap.results, expected)
