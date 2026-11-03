from unittest import TestCase
from copy import deepcopy
import gzip
import tempfile
import os
import json

from pydantic import ValidationError

from pisa_utils.post_process_parsers import (
    PostProcessComplexTable,
    PostProcessComponentsTable,
    PostProcessInterfaceDetailsList,
    PostProcessInterfaceTable,
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

        self.path_interface_jsons = os.path.join(
            self.base_data_path, "mock_data/post_process_interface_table/interfaces"
        )
        self.interface_jsons = {}
        for interface_id in (1, 2, 3, 4):
            with open(
                os.path.join(
                    self.path_interface_jsons, f"interface_{interface_id}.json"
                ),
                "r",
            ) as f:
                self.interface_jsons[interface_id] = json.load(f)

    def _write_interface_jsons(self, directory, interfaces, compressed=False):
        """Write individual interface fixtures in the requested compression format."""
        os.makedirs(directory)
        extension = ".json.gz" if compressed else ".json"
        open_json = gzip.open if compressed else open
        for interface_id, data in interfaces.items():
            with open_json(
                os.path.join(directory, f"interface_{interface_id}" + extension), "wt"
            ) as f:
                json.dump(data, f)

    def test_parse_interface_table(self):
        """
        Read a summary and individual interfaces from JSON or gzip JSON files.

        Extract all three bond counts, replacing stale summary counts, and
        compare the complete output including rounded area, energy, and scores.
        """
        expected_data = deepcopy(self.data)
        expected_data["interface_types"][0]["interfaces"][0].update(
            int_area=1617.4,
            int_solv_energy=-23.1,
            pvalue=0.661,
            css=0.382,
            n_h_bonds=5,
            n_salt_bridges=2,
            n_ss_bonds=1,
        )

        for compressed in (False, True):
            with (
                self.subTest(compressed=compressed),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                extension = ".json.gz" if compressed else ".json"
                input_json_path = os.path.join(tmpdir, "interface_summary" + extension)
                output_json_path = os.path.join(tmpdir, "output" + extension)
                path_interface_jsons = os.path.join(tmpdir, "interfaces")
                self._write_interface_jsons(
                    path_interface_jsons, self.interface_jsons, compressed=compressed
                )
                open_json = gzip.open if compressed else open
                with open_json(input_json_path, "wt") as f:
                    json.dump(self.data, f)

                parser = PostProcessInterfaceTable(
                    input_summary_interfaces_path=input_json_path,
                    path_interface_jsons=path_interface_jsons,
                    output_json_path=output_json_path,
                    compressed=compressed,
                )
                parser.parse()

                with open_json(output_json_path, "rt") as f:
                    output_data = json.load(f)
                self.assertDictEqual(output_data, expected_data)

    def test_parse_interface_table_bond_counts(self):
        """
        Match bond counts to interface IDs and default missing counts to zero.

        Cover distinct counts in two files, omitted bond sections, and sections
        lacking n_bonds, using summaries that do not already have bond counts.
        """
        template = deepcopy(self.data["interface_types"][0]["interfaces"][0])
        for field in ("n_h_bonds", "n_salt_bridges", "n_ss_bonds"):
            template.pop(field)
        data = {
            "n_interfaces": 4,
            "interface_types": [
                {
                    "int_type": 1,
                    "interfaces": [
                        {**template, "interface_id": interface_id}
                        for interface_id in (2, 4, 1, 3)
                    ],
                }
            ],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            input_json_path = os.path.join(tmpdir, "interface_summary.json")
            output_json_path = os.path.join(tmpdir, "output.json")
            with open(input_json_path, "w") as f:
                json.dump(data, f)

            parser = PostProcessInterfaceTable(
                input_summary_interfaces_path=input_json_path,
                path_interface_jsons=self.path_interface_jsons,
                output_json_path=output_json_path,
            )
            parser.parse()

            with open(output_json_path, "r") as f:
                output_data = json.load(f)

            self.assertEqual(output_data["n_interfaces"], 4)
            self.assertDictEqual(
                {
                    row["interface_id"]: (
                        row["n_h_bonds"],
                        row["n_salt_bridges"],
                        row["n_ss_bonds"],
                    )
                    for row in output_data["interface_types"][0]["interfaces"]
                },
                {1: (5, 2, 1), 2: (11, 4, 2), 3: (0, 0, 0), 4: (0, 0, 0)},
            )

    def test_parse_interface_table_ordering(self):
        """
        Order groups by int_type and interfaces within each group by interface_id.

        Start with shuffled groups and rows, including IDs above nine, to check
        ascending numeric ordering and keep bond counts with the correct rows.
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
            path_interface_jsons = os.path.join(tmpdir, "interfaces")
            self._write_interface_jsons(
                path_interface_jsons,
                {
                    row["interface_id"]: {
                        "interface": {"h_bonds": {"n_bonds": row["interface_id"]}}
                    }
                    for group in data["interface_types"]
                    for row in group["interfaces"]
                },
            )
            with open(input_json_path, "w") as f:
                json.dump(data, f)

            parser = PostProcessInterfaceTable(
                input_summary_interfaces_path=input_json_path,
                path_interface_jsons=path_interface_jsons,
                output_json_path=output_json_path,
            )
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
            for group in output_data["interface_types"]:
                for row in group["interfaces"]:
                    self.assertEqual(row["n_h_bonds"], row["interface_id"])

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

                parser = PostProcessInterfaceTable(
                    input_summary_interfaces_path=input_json_path,
                    path_interface_jsons=self.path_interface_jsons,
                    output_json_path=output_json_path,
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

                parser = PostProcessInterfaceTable(
                    input_summary_interfaces_path=input_json_path,
                    path_interface_jsons=self.path_interface_jsons,
                    output_json_path=output_json_path,
                )
                with self.assertRaises(expected_error):
                    parser.parse()

                self.assertFalse(os.path.exists(output_json_path))

    def test_parse_interface_table_missing_interface_json(self):
        """Raise for a missing individual interface file without writing output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_json_path = os.path.join(tmpdir, "interface_summary.json")
            output_json_path = os.path.join(tmpdir, "output.json")
            with open(input_json_path, "w") as f:
                json.dump(self.data, f)

            parser = PostProcessInterfaceTable(
                input_summary_interfaces_path=input_json_path,
                path_interface_jsons=tmpdir,
                output_json_path=output_json_path,
            )
            with self.assertRaises(FileNotFoundError) as context:
                parser.parse()

            self.assertEqual(
                context.exception.filename, os.path.join(tmpdir, "interface_1.json")
            )
            self.assertFalse(os.path.exists(output_json_path))


class TestPostProcessComponentsTable(PostProcessTestCase):
    def setUp(self):
        super().setUp()
        input_json_path = os.path.join(
            self.base_data_path,
            "mock_data/post_process_components_table/monomers_extended.json",
        )
        expected_json_path = os.path.join(
            self.base_data_path,
            "mock_data/model_inputs/post_processor_models/components_table.json",
        )
        with open(input_json_path, "r") as f:
            self.data = json.load(f)
        with open(expected_json_path, "r") as f:
            self.expected_data = json.load(f)

    def test_parse_components_table(self):
        """
        Group shuffled monomers and write the complete table as JSON or gzip JSON.

        Check numeric type ordering, component ID ordering within shared types,
        rounding, ligand fields, and removal of type IDs from individual rows.
        """
        for compressed in (False, True):
            with (
                self.subTest(compressed=compressed),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                extension = ".json.gz" if compressed else ".json"
                input_json_path = os.path.join(tmpdir, "monomers_extended" + extension)
                output_json_path = os.path.join(tmpdir, "components_table" + extension)
                open_json = gzip.open if compressed else open
                with open_json(input_json_path, "wt") as f:
                    json.dump(self.data, f)

                parser = PostProcessComponentsTable(
                    input_json_path, output_json_path, compressed=compressed
                )
                parser.parse()

                with open_json(output_json_path, "rt") as f:
                    output_data = json.load(f)
                self.assertListEqual(output_data, self.expected_data)

    def test_parse_components_table_empty(self):
        """Write an empty JSON list when the monomer input contains no components."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_json_path = os.path.join(tmpdir, "monomers_extended.json")
            output_json_path = os.path.join(tmpdir, "components_table.json")
            with open(input_json_path, "w") as f:
                json.dump({"components": []}, f)

            parser = PostProcessComponentsTable(input_json_path, output_json_path)
            parser.parse()

            with open(output_json_path, "r") as f:
                self.assertEqual(json.load(f), [])

    def test_parse_components_table_invalid_input(self):
        """Reject missing files, malformed JSON, and invalid monomers without output."""
        missing_type = deepcopy(self.data)
        missing_type["components"][0].pop("molecule_type_id")
        invalid_area = deepcopy(self.data)
        invalid_area["components"][0]["asa"] = "invalid"
        unexpected_field = deepcopy(self.data)
        unexpected_field["components"][0]["unexpected"] = True
        cases = (
            ("missing_file", None, FileNotFoundError),
            ("malformed_json", "{", json.JSONDecodeError),
            ("missing_components", "{}", ValidationError),
            ("null_components", json.dumps({"components": None}), ValidationError),
            ("missing_type", json.dumps(missing_type), ValidationError),
            ("invalid_area", json.dumps(invalid_area), ValidationError),
            ("unexpected_field", json.dumps(unexpected_field), ValidationError),
        )
        for case, contents, expected_error in cases:
            with (
                self.subTest(case=case),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                input_json_path = os.path.join(tmpdir, "monomers_extended.json")
                output_json_path = os.path.join(tmpdir, "components_table.json")
                if contents is not None:
                    with open(input_json_path, "w") as f:
                        f.write(contents)

                parser = PostProcessComponentsTable(input_json_path, output_json_path)
                with self.assertRaises(expected_error):
                    parser.parse()

                self.assertFalse(os.path.exists(output_json_path))
