from unittest import TestCase
from copy import deepcopy
import gzip
import tempfile
import os
import json

from pydantic import ValidationError

from pisa_utils.post_process_parsers import (
    PostProcessComplexTable,
    PostProcessInterfaceDetailsList,
    PostProcessInterfaceTableParser,
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


class TestPostProcessInterfaceTableParser(PostProcessTestCase):
    def setUp(self):
        super().setUp()
        input_json_path = os.path.join(
            self.base_data_path,
            "mock_data/model_inputs/post_processor_models/interface_table.json",
        )
        with open(input_json_path, "r") as f:
            self.data = json.load(f)

    def test_parse_interface_table(self):
        """
        Read a complete summary and write rounded values to JSON or gzip JSON.

        Compare the full output to ensure all fields and groups are preserved
        while area, solvation energy, pvalue, and CSS are rounded correctly.
        """
        expected_data = deepcopy(self.data)
        expected_data["interface_types"][0]["interfaces"][0].update(
            int_area=1617.4,
            int_solv_energy=-23.1,
            pvalue=0.661,
            css=0.382,
        )

        for compressed in (False, True):
            with (
                self.subTest(compressed=compressed),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                extension = ".json.gz" if compressed else ".json"
                input_json_path = os.path.join(tmpdir, "interface_summary" + extension)
                output_json_path = os.path.join(tmpdir, "output" + extension)
                open_json = gzip.open if compressed else open
                with open_json(input_json_path, "wt") as f:
                    json.dump(self.data, f)

                parser = PostProcessInterfaceTableParser(
                    input_json_path, output_json_path, compressed=compressed
                )
                parser.parse()

                with open_json(output_json_path, "rt") as f:
                    output_data = json.load(f)
                self.assertDictEqual(output_data, expected_data)

    def test_parse_interface_table_ordering(self):
        """
        Order groups by int_type and interfaces within each group by interface_id.

        Start with shuffled groups and rows, including IDs above nine, to check
        ascending numeric ordering at both levels of the saved JSON.
        """
        template = self.data["interface_types"][0]["interfaces"][0]
        data = {"n_interfaces": 9, "interface_types": []}
        for int_type, interface_ids in (
            (2, [12, 3, 7]),
            (10, [11, 5, 8]),
            (1, [10, 2, 6]),
        ):
            data["interface_types"].append(
                {
                    "int_type": int_type,
                    "interfaces": [
                        {**template, "interface_id": interface_id}
                        for interface_id in interface_ids
                    ],
                }
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            input_json_path = os.path.join(tmpdir, "interface_summary.json")
            output_json_path = os.path.join(tmpdir, "output.json")
            with open(input_json_path, "w") as f:
                json.dump(data, f)

            parser = PostProcessInterfaceTableParser(input_json_path, output_json_path)
            parser.parse()

            with open(output_json_path, "r") as f:
                output_data = json.load(f)

            self.assertEqual(output_data["n_interfaces"], 9)
            self.assertListEqual(
                [
                    (
                        group["int_type"],
                        [row["interface_id"] for row in group["interfaces"]],
                    )
                    for group in output_data["interface_types"]
                ],
                [(1, [2, 6, 10]), (2, [3, 7, 12]), (10, [5, 8, 11])],
            )

    def test_parse_interface_table_empty(self):
        """
        Write valid summaries containing no interfaces without adding rows.

        Cover both an empty interface_types list and a group with an empty
        interfaces list, preserving the input structure and zero count.
        """
        for groups in ([], [{"int_type": 1, "interfaces": []}]):
            with self.subTest(groups=groups), tempfile.TemporaryDirectory() as tmpdir:
                data = {"n_interfaces": 0, "interface_types": groups}
                input_json_path = os.path.join(tmpdir, "interface_summary.json")
                output_json_path = os.path.join(tmpdir, "output.json")
                with open(input_json_path, "w") as f:
                    json.dump(data, f)

                parser = PostProcessInterfaceTableParser(
                    input_json_path, output_json_path
                )
                parser.parse()

                with open(output_json_path, "r") as f:
                    output_data = json.load(f)
                self.assertDictEqual(output_data, data)

    def test_parse_interface_table_invalid_input(self):
        """
        Propagate input errors without creating an output file.

        Cover a missing input file, malformed JSON, and a summary missing the
        required n_interfaces field, checking the corresponding exception type.
        """
        cases = (
            (None, FileNotFoundError),
            ("{", json.JSONDecodeError),
            (json.dumps({"interface_types": []}), ValidationError),
        )
        for contents, expected_error in cases:
            with (
                self.subTest(expected_error=expected_error),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                input_json_path = os.path.join(tmpdir, "interface_summary.json")
                output_json_path = os.path.join(tmpdir, "output.json")
                if contents is not None:
                    with open(input_json_path, "w") as f:
                        f.write(contents)

                parser = PostProcessInterfaceTableParser(
                    input_json_path, output_json_path
                )
                with self.assertRaises(expected_error):
                    parser.parse()

                self.assertFalse(os.path.exists(output_json_path))
