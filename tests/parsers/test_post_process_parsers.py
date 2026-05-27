from unittest import TestCase
import tempfile
import os
import json

from pisa_utils.post_process_parsers import (
    PostProcessComplexTable,
    PostProcessInterfaceDetailsList,
)


class PostProcessTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.maxDiff = None
        self.base_data_path = "tests/data/"


class TestPostProcessParsers(PostProcessTestCase):
    def test_parse_complex_table(self):
        """
        Test a parse with a set of real, truncated data from 3hax.
        """

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


class TestPostProcessInterfaceDetailsList(PostProcessTestCase):
    def test_parse_interface_details_list(self):
        """
        Test a parse with a set of real, truncated data from 3hax.
        """

        expected_output_path = os.path.join(
            self.base_data_path,
            "expected_output",
            "post_processed_jsons",
            "3hax_multi_interface_details_list.json",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path_interfaces = os.path.join(
                self.base_data_path,
                "expected_output",
                "interfaces",
                "3hax_interfaces_multi",
            )
            path_monomers_json = os.path.join(
                self.base_data_path,
                "mock_data/list_results/3hax/monomers_extended.json",
            )

            output_json_path = os.path.join(tmpdir, "output.json")

            # Run parser
            parser = PostProcessInterfaceDetailsList(
                path_interfaces=path_interfaces,
                path_monomers_json=path_monomers_json,
                output_json_path=output_json_path,
            )
            parser.parse()

            # Checks
            with open(expected_output_path, "r") as f:
                expected_data = json.load(f)
            with open(output_json_path, "r") as f:
                output_data = json.load(f)

            self.assertListEqual(expected_data, output_data)
