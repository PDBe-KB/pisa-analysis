import json
import os
import unittest
from unittest import TestCase

from pisa_utils import analyze

test_data_path = os.path.join("tests", "data")


class TestJsonFile(TestCase):
    def test_json_output_equal_result(self):

        ref_data_file = os.path.join(test_data_path, "test_result.json")
        with open(ref_data_file) as test_file:
            ref_json = json.load(test_file)

        ap = analyze.AnalysePisa(
            pdb_id="6nxr",
            assembly_id="1",
            output_dir=test_data_path,
            result_json_file="output.json",
            input_dir=test_data_path,
            input_updated_cif=None,
            input_cif_file=None,
        )

        ap.process_pisa_xml()
        ap.set_results()
        ap.save_to_json()

        output_data_file = os.path.join(test_data_path, "output.json")
        with open(output_data_file) as output_file:
            out_json = json.load(output_file)

        self.assertEqual(ref_json, out_json)


if __name__ == "__main__":
    unittest.main()
