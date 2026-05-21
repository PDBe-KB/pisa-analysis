from unittest import TestCase
import tempfile
import os
import json

from pisa_utils.post_process_parsers import PostProcessComplexTable


class TestPostProcessParsers(TestCase):
    def setUp(self):
        super().setUp()
        self.maxDiff = None
        self.base_data_path = "tests/data/"

    def test_parse_complex_table(self):

        expected_output_path = os.path.join(
            self.base_data_path,
            "expected_output",
            "post_processed_jsons",
            "3hax_assembly_multi_asmset_post_proc.json",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            input_json_path = os.path.join(
                self.base_data_path,
                "expected_output",
                "3hax_assembly_multi_asmset.json",
            )
            output_json_path = os.path.join(tmpdir, "output.json")

            # Run parser
            parser = PostProcessComplexTable(input_json_path, output_json_path)
            parser.parse()

            # Checks
            with open(expected_output_path, "r") as f:
                expected_data = json.load(f)
            with open(output_json_path, "r") as f:
                output_data = json.load(f)

            self.assertListEqual(expected_data, output_data)
