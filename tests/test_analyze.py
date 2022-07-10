import os
import xml.etree.ElementTree as ET
from unittest import TestCase
from unittest.mock import patch

from pisa_utils.analyze import AnalysePisa


class TestAnalyzePisa(TestCase):
    @patch("pisa_utils.analyze.get_bond_dict")
    @patch(
        "pisa_utils.analyze.parse_xml_file",
        side_effect=[
            ET.parse(
                os.path.join(".", "tests", "data", "mocks", "assembly.xml")
            ).getroot(),
            ET.parse(
                os.path.join(".", "tests", "data", "mocks", "interfaces.xml")
            ).getroot(),
        ],
    )
    def test_process_pisa_xml(self, bonds, xmls):
        bonds.return_value = {}
        ap = AnalysePisa(
            pdb_id="6nxr",
            assembly_id="1",
            output_dir=os.path.join("tests", "data", "mocks"),
            result_json_file=".",
            input_dir=".",
            input_updated_cif=".",
            input_cif_file=".",
        )
        expected = {
            "assembly_status": "Ok",
            "status": "Ok",
            "num_interfaces": "5",
            "non_ligand_interface_count": 0,
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

        self.assertEqual(ap.process_pisa_xml(), expected)

    # def test_set_results(self):
    #     ap = AnalysePisa(
    #         pdb_id='6nxr',
    #         assembly_id='1',
    #         output_dir='.',
    #         result_json_file='.',
    #         input_dir='.',
    #         input_updated_cif='.',
    #         input_cif_file='.',
    #     )
    #
    #     ap.set_results()
    #     self.assertEqual(ap.results['PISA'], {})
