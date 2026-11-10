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
    PostProcessInterfaceTable,
    PostProcessInterfaceComponentResidueTable,
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

        input_json_path = os.path.join(
            self.base_data_path,
            "expected_output",
            "3hax_assembly_multi_asmset.json",
        )
        with open(expected_output_path) as f:
            expected_data = json.load(f)

        for compressed in (False, True):
            with (
                self.subTest(compressed=compressed),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                source_path = input_json_path
                output_json_path = os.path.join(tmpdir, "output.json")
                open_json = gzip.open if compressed else open
                if compressed:
                    source_path = os.path.join(tmpdir, "assemblies.json.gz")
                    output_json_path += ".gz"
                    with open(input_json_path) as source:
                        with gzip.open(source_path, "wt") as target:
                            target.write(source.read())

                parser = PostProcessComplexTable(
                    source_path, output_json_path, compressed=compressed
                )
                parser.parse()

                with open_json(output_json_path, "rt") as f:
                    output_data = json.load(f)

                self.assertListEqual(expected_data, output_data)

    def test_parse_complex_table_with_missing_groups(self):
        """Handle ASU-only, PQS-only, and empty results without nested groups."""
        with open(
            os.path.join(
                self.base_data_path, "expected_output/3hax_assembly_multi_asmset.json"
            )
        ) as f:
            assemblies = json.load(f)
        with open(
            os.path.join(
                self.base_data_path,
                "expected_output/post_processed_jsons/3hax_assembly_multi_asmset_post_proc.json",
            )
        ) as f:
            expected = json.load(f)

        cases = (
            ({"asu_complex": assemblies["asu_complex"]}, expected[:1]),
            ({"pqs_sets": assemblies["pqs_sets"]}, expected[1:]),
            ({"pqs_sets": [{"pqs_set_id": 1, "complexes": []}]}, []),
            ({}, []),
        )
        for assemblies, expected in cases:
            with (
                self.subTest(assemblies=assemblies),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                input_path = os.path.join(tmpdir, "assemblies.json")
                output_path = os.path.join(tmpdir, "complex_table.json")
                with open(input_path, "w") as f:
                    json.dump(assemblies, f)

                PostProcessComplexTable(input_path, output_path).parse()

                with open(output_path) as f:
                    self.assertEqual(json.load(f), expected)


class TestPostProcessInterfaceDetailsList(PostProcessTestCase):
    def setUp(self):
        super().setUp()
        self.path_expected = os.path.join(
            self.base_data_path,
            "expected_output/post_processed_jsons/3hax_multi_interface_details",
        )
        self.path_interfaces = os.path.join(
            self.base_data_path, "expected_output/interfaces/3hax_interfaces_multi"
        )
        self.path_monomers = os.path.join(
            self.base_data_path, "mock_data/list_results/3hax/monomers_extended.json"
        )

    def test_parse_interface_details_files(self):
        """Join real data into individual JSON or gzip files, creating the folder."""
        for compressed in (False, True):
            with (
                self.subTest(compressed=compressed),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                path_interfaces = os.path.join(tmpdir, "interfaces")
                output_json_dir = os.path.join(tmpdir, "details", "individual")
                os.makedirs(path_interfaces)
                extension = ".json.gz" if compressed else ".json"
                open_json = gzip.open if compressed else open
                for interface_id in (3, 1, 2):
                    filename = f"interface_{interface_id}.json"
                    with open(os.path.join(self.path_interfaces, filename)) as source:
                        with open_json(
                            os.path.join(
                                path_interfaces, f"interface_{interface_id}" + extension
                            ),
                            "wt",
                        ) as target:
                            target.write(source.read())
                # Ignore unrelated files and directories resembling interface files.
                for name in ("interface_notes", "interface_1_details"):
                    with open(
                        os.path.join(path_interfaces, name + extension), "w"
                    ) as f:
                        f.write("not JSON")
                os.mkdir(os.path.join(path_interfaces, "interface_4" + extension))

                path_monomers = self.path_monomers
                if compressed:
                    path_monomers = os.path.join(tmpdir, "monomers.json.gz")
                    with open(self.path_monomers) as source:
                        with gzip.open(path_monomers, "wt") as target:
                            target.write(source.read())

                parser = PostProcessInterfaceDetailsList(
                    path_interfaces=path_interfaces,
                    path_monomers_json=path_monomers,
                    output_json_dir=output_json_dir,
                    compressed=compressed,
                )
                parser.parse()

                self.assertEqual(
                    set(os.listdir(output_json_dir)),
                    {
                        f"interface_{interface_id}" + extension
                        for interface_id in (1, 2, 3)
                    },
                )
                for interface_id in (1, 2, 3):
                    with open(
                        os.path.join(
                            self.path_expected, f"interface_{interface_id}.json"
                        )
                    ) as f:
                        expected = json.load(f)
                    with open_json(
                        os.path.join(
                            output_json_dir, f"interface_{interface_id}" + extension
                        ),
                        "rt",
                    ) as f:
                        self.assertDictEqual(json.load(f), expected)

    def test_missing_or_empty_directory(self):
        """Skip absent or empty interface directories without creating output."""
        for missing in (False, True):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as tmpdir:
                path_interfaces = os.path.join(tmpdir, "interfaces")
                output_json_dir = os.path.join(tmpdir, "details")
                if not missing:
                    os.makedirs(path_interfaces)

                parser = PostProcessInterfaceDetailsList(
                    path_interfaces=path_interfaces,
                    path_monomers_json=self.path_monomers,
                    output_json_dir=output_json_dir,
                )
                with self.assertLogs(
                    "pisa_utils.post_process_parsers", level="WARNING"
                ):
                    parser.parse()
                self.assertFalse(os.path.exists(output_json_dir))

    def test_invalid_monomer_data(self):
        """Reject missing or mismatched monomers before writing the interface file."""
        with open(os.path.join(self.path_interfaces, "interface_1.json")) as f:
            interface = json.load(f)
        component_id = interface["interface"]["molecules"][0]["component_id"]

        for invalid in ("missing", "mismatched_chain"):
            with self.subTest(invalid=invalid), tempfile.TemporaryDirectory() as tmpdir:
                with open(self.path_monomers) as f:
                    monomers = json.load(f)
                if invalid == "missing":
                    monomers["components"] = [
                        component
                        for component in monomers["components"]
                        if component["component_id"] != component_id
                    ]
                else:
                    next(
                        component
                        for component in monomers["components"]
                        if component["component_id"] == component_id
                    )["auth_asym_id"] = "other"

                path_interfaces = os.path.join(tmpdir, "interfaces")
                os.makedirs(path_interfaces)
                with open(os.path.join(path_interfaces, "interface_1.json"), "w") as f:
                    json.dump(interface, f)
                path_monomers = os.path.join(tmpdir, "monomers.json")
                with open(path_monomers, "w") as f:
                    json.dump(monomers, f)
                output_json_dir = os.path.join(tmpdir, "details")
                parser = PostProcessInterfaceDetailsList(
                    path_interfaces=path_interfaces,
                    path_monomers_json=path_monomers,
                    output_json_dir=output_json_dir,
                )
                with self.assertRaises(ValueError):
                    parser.parse()
                self.assertFalse(
                    os.path.exists(os.path.join(output_json_dir, "interface_1.json"))
                )


class TestPostProcessInterfaceComponentResidueTable(PostProcessTestCase):
    def setUp(self):
        super().setUp()
        self.data = {
            "interface_id": 2,
            "interface": {
                "molecules": [
                    {
                        "auth_asym_id": "A",
                        "residues": {
                            "residues": [
                                {
                                    "residue_serial_number": 1,
                                    "auth_comp_id": "ARG",
                                    "auth_seq_id": 11,
                                    "label_seq_id": 8,
                                    "ins_code": None,
                                    "bonds": "HS",
                                    "asa": 158.712,
                                    "bsa": 123.567,
                                    "solv_energy": -1.23456,
                                },
                                {
                                    "auth_comp_id": "GLY",
                                    "auth_seq_id": -2,
                                    "asa": 100.0,
                                    "bsa": 0.0,
                                    "solv_energy": 0.0,
                                },
                            ]
                        },
                    },
                    {
                        "auth_asym_id": "A",
                        "residues": {
                            "residues": [
                                {
                                    "auth_comp_id": "CYS",
                                    "auth_seq_id": 11,
                                    "bonds": None,
                                    "asa": 0.0,
                                    "bsa": 0.0,
                                    "solv_energy": 0.0,
                                }
                            ]
                        },
                    },
                ]
            },
        }

    def test_parse_residue_tables(self):
        """Keep both sides, order, and all accessibility classes in JSON or gzip."""
        expected = {
            "structure_1_residues": [
                {
                    "auth_asym_id": "A",
                    "label_comp_id": "ARG",
                    "auth_seq_id": 11,
                    "bonds": "HS",
                    "asa": 158.7,
                    "bsa": 123.6,
                    "pc_residue_buried": 78.0,
                    "solv_energy": -1.235,
                    "is_interfacing_residue": True,
                    "is_solvent_accessible_residue": False,
                    "is_inaccessible_residue": False,
                },
                {
                    "auth_asym_id": "A",
                    "label_comp_id": "GLY",
                    "auth_seq_id": -2,
                    "bonds": None,
                    "asa": 100.0,
                    "bsa": 0.0,
                    "pc_residue_buried": 0.0,
                    "solv_energy": 0.0,
                    "is_interfacing_residue": False,
                    "is_solvent_accessible_residue": True,
                    "is_inaccessible_residue": False,
                },
            ],
            "structure_2_residues": [
                {
                    "auth_asym_id": "A",
                    "label_comp_id": "CYS",
                    "auth_seq_id": 11,
                    "bonds": None,
                    "asa": 0.0,
                    "bsa": 0.0,
                    "pc_residue_buried": 0.0,
                    "solv_energy": 0.0,
                    "is_interfacing_residue": False,
                    "is_solvent_accessible_residue": False,
                    "is_inaccessible_residue": True,
                }
            ],
        }

        for compressed in (False, True):
            with (
                self.subTest(compressed=compressed),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                path_interfaces = os.path.join(tmpdir, "interfaces")
                output_json_dir = os.path.join(tmpdir, "tables")
                os.makedirs(path_interfaces)
                extension = ".json.gz" if compressed else ".json"
                open_json = gzip.open if compressed else open
                empty_interface = {
                    "interface_id": 10,
                    "interface": {
                        "molecules": [
                            {"auth_asym_id": "B", "residues": {"residues": []}},
                            {"auth_asym_id": "C", "residues": {"residues": []}},
                        ]
                    },
                }
                for data in (self.data, empty_interface):
                    with open_json(
                        os.path.join(
                            path_interfaces,
                            f"interface_{data['interface_id']}" + extension,
                        ),
                        "wt",
                    ) as f:
                        json.dump(data, f)
                # Ignore unrelated files and directories resembling interface files.
                for name in ("summary", "interface_notes", "interface_2_residue_table"):
                    with open(
                        os.path.join(path_interfaces, name + extension), "w"
                    ) as f:
                        f.write("not JSON")
                os.mkdir(os.path.join(path_interfaces, "interface_3" + extension))

                parser = PostProcessInterfaceComponentResidueTable(
                    path_interfaces=path_interfaces,
                    output_json_dir=output_json_dir,
                    compressed=compressed,
                )
                parser.parse()

                self.assertEqual(
                    set(os.listdir(output_json_dir)),
                    {
                        "interface_2_residue_table" + extension,
                        "interface_10_residue_table" + extension,
                    },
                )
                for interface_id, expected_data in (
                    (2, expected),
                    (10, {"structure_1_residues": [], "structure_2_residues": []}),
                ):
                    with open_json(
                        os.path.join(
                            output_json_dir,
                            f"interface_{interface_id}_residue_table" + extension,
                        ),
                        "rt",
                    ) as f:
                        self.assertDictEqual(json.load(f), expected_data)

    def test_parse_residue_tables_rejects_removed_energy_alias(self):
        """Require solv_energy in input residues without writing invalid tables."""
        residue = self.data["interface"]["molecules"][0]["residues"]["residues"][0]
        residue["solv_en"] = residue.pop("solv_energy")

        with tempfile.TemporaryDirectory() as tmpdir:
            with open(os.path.join(tmpdir, "interface_2.json"), "w") as f:
                json.dump(self.data, f)

            parser = PostProcessInterfaceComponentResidueTable(tmpdir, tmpdir)
            with self.assertRaises(KeyError) as context:
                parser.parse()

            self.assertEqual(context.exception.args, ("solv_energy",))
            self.assertFalse(
                os.path.exists(os.path.join(tmpdir, "interface_2_residue_table.json"))
            )

    def test_parse_real_interfaces(self):
        """Read existing interface JSONs, including a single-residue ligand."""
        path_interfaces = os.path.join(
            self.base_data_path, "expected_output/interfaces/af3_heterodimer_with_atp"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            parser = PostProcessInterfaceComponentResidueTable(path_interfaces, tmpdir)
            parser.parse()

            self.assertEqual(len(os.listdir(tmpdir)), 3)
            with open(os.path.join(tmpdir, "interface_3_residue_table.json")) as f:
                table = json.load(f)
            self.assertEqual(len(table["structure_1_residues"]), 30)
            self.assertEqual(len(table["structure_2_residues"]), 1)
            ligand = table["structure_2_residues"][0]
            self.assertEqual(ligand["auth_asym_id"], "C")
            self.assertEqual(ligand["label_comp_id"], "ATP")
            self.assertEqual(ligand["bonds"], "H")
            self.assertEqual(ligand["pc_residue_buried"], 14.0)
            self.assertTrue(ligand["is_interfacing_residue"])

    def test_missing_or_empty_directory(self):
        """Match other directory parsers by skipping absent or empty inputs."""
        for missing in (False, True):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as tmpdir:
                path_interfaces = os.path.join(tmpdir, "interfaces")
                output_json_dir = os.path.join(tmpdir, "tables")
                if not missing:
                    os.makedirs(path_interfaces)

                parser = PostProcessInterfaceComponentResidueTable(
                    path_interfaces, output_json_dir
                )
                with self.assertLogs(
                    "pisa_utils.post_process_parsers", level="WARNING"
                ):
                    parser.parse()
                self.assertFalse(os.path.exists(output_json_dir))

    def test_invalid_interface(self):
        """Reject incorrect component counts and invalid residue metrics."""
        for molecule_count in (0, 1, 3):
            with (
                self.subTest(molecule_count=molecule_count),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                data = deepcopy(self.data)
                data["interface"]["molecules"] = [
                    data["interface"]["molecules"][0]
                ] * molecule_count
                with open(os.path.join(tmpdir, "interface_2.json"), "w") as f:
                    json.dump(data, f)

                parser = PostProcessInterfaceComponentResidueTable(tmpdir, tmpdir)
                with self.assertRaisesRegex(ValueError, "exactly two molecules"):
                    parser.parse()
                self.assertFalse(
                    os.path.exists(
                        os.path.join(tmpdir, "interface_2_residue_table.json")
                    )
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            self.data["interface"]["molecules"][0]["residues"]["residues"][0]["asa"] = (
                "bad"
            )
            with open(os.path.join(tmpdir, "interface_2.json"), "w") as f:
                json.dump(self.data, f)

            parser = PostProcessInterfaceComponentResidueTable(tmpdir, tmpdir)
            with self.assertRaises(ValidationError):
                parser.parse()
            self.assertFalse(
                os.path.exists(os.path.join(tmpdir, "interface_2_residue_table.json"))
            )


class TestPostProcessInterfaceTableParser(PostProcessTestCase):
    def setUp(self):
        super().setUp()
        input_json_path = os.path.join(
            self.base_data_path,
            "mock_data/post_process_interface_table/interface_summary.json",
        )
        with open(input_json_path, "r") as f:
            self.data = json.load(f)

        expected_json_path = os.path.join(
            self.base_data_path,
            "mock_data/model_inputs/post_processor_models/interface_table.json",
        )
        with open(expected_json_path) as f:
            self.expected_data = json.load(f)

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
        expected_data = deepcopy(self.expected_data)
        expected_data[0].update(
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
                self.assertListEqual(output_data, expected_data)

    def test_parse_interface_table_rejects_removed_aliases(self):
        """Reject summary rows using removed aliases without writing output."""
        for field, alias in (
            ("auth_asym_id_2", "chain_id"),
            ("int_solv_energy", "int_solv_en"),
        ):
            with self.subTest(alias=alias), tempfile.TemporaryDirectory() as tmpdir:
                data = deepcopy(self.data)
                row = data["interface_types"][0]["interfaces"][0]
                row[alias] = row.pop(field)
                input_json_path = os.path.join(tmpdir, "interface_summary.json")
                output_json_path = os.path.join(tmpdir, "output.json")
                with open(input_json_path, "w") as f:
                    json.dump(data, f)

                parser = PostProcessInterfaceTable(
                    input_summary_interfaces_path=input_json_path,
                    path_interface_jsons=self.path_interface_jsons,
                    output_json_path=output_json_path,
                )
                with self.assertRaises(ValidationError) as context:
                    parser.parse()

                self.assertEqual(
                    {
                        (error["loc"], error["type"])
                        for error in context.exception.errors()
                    },
                    {((0, field), "missing"), ((0, alias), "extra_forbidden")},
                )
                self.assertFalse(os.path.exists(output_json_path))

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

            self.assertEqual(len(output_data), 4)
            self.assertDictEqual(
                {
                    row["interface_id"]: (
                        row["n_h_bonds"],
                        row["n_salt_bridges"],
                        row["n_ss_bonds"],
                    )
                    for row in output_data
                },
                {1: (5, 2, 1), 2: (11, 4, 2), 3: (0, 0, 0), 4: (0, 0, 0)},
            )

    def test_parse_interface_table_ordering(self):
        """
        Order flat rows by int_type, then interface_id.

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

            self.assertEqual(len(output_data), 9)
            self.assertListEqual(
                [(row["int_type"], row["interface_id"]) for row in output_data],
                [
                    (1, 2),
                    (1, 6),
                    (1, 10),
                    (2, 3),
                    (2, 7),
                    (2, 12),
                    (10, 5),
                    (10, 8),
                    (10, 11),
                ],
            )
            for row in output_data:
                self.assertEqual(row["n_h_bonds"], row["interface_id"])

    def test_parse_interface_table_empty(self):
        """
        Write valid summaries containing no interfaces without adding rows.

        Cover both an empty interface_types list and a group with an empty
        interfaces list, writing an empty JSON list in both cases.
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
                self.assertListEqual(output_data, [])

    def test_parse_interface_table_invalid_input(self):
        """
        Propagate input errors without creating an output file.

        Cover a missing input file, malformed JSON, and a summary row missing a
        required field, checking the corresponding exception type.
        """
        invalid_summary = deepcopy(self.data)
        del invalid_summary["interface_types"][0]["interfaces"][0]["auth_asym_id_1"]
        cases = (
            (None, FileNotFoundError),
            ("{", json.JSONDecodeError),
            (json.dumps(invalid_summary), ValidationError),
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
