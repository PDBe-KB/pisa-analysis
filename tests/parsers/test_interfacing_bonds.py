from copy import deepcopy
import gzip
import json
import os
import tempfile
from unittest import TestCase

from pydantic import ValidationError

from pisa_utils.post_process_parsers import PostProcessInterfacingBondsTable


class TestPostProcessInterfacingBondsTable(TestCase):
    def setUp(self):
        self.maxDiff = None
        self.bond = {
            "auth_asym_id_1": "B",
            "auth_comp_id_1": "ARG",
            "auth_seq_id_1": -2,
            "auth_atom_id_1": "NH1",
            "auth_asym_id_2": "A",
            "auth_comp_id_2": "GLU",
            "auth_seq_id_2": 105,
            "auth_atom_id_2": "OE2",
            "dist": 2.87654,
            # Source fields outside the output schema must be omitted.
            "label_asym_id_1": "C",
            "label_seq_id_1": 1,
            "orig_label_asym_id_2": None,
            "inscode_1": "A",
        }
        self.expected_bond = {
            "auth_asym_id_1": "B",
            "label_comp_id_1": "ARG",
            "auth_seq_id_1": -2,
            "auth_atom_id_1": "NH1",
            "auth_asym_id_2": "A",
            "label_comp_id_2": "GLU",
            "auth_seq_id_2": 105,
            "auth_atom_id_2": "OE2",
            "dist": 2.877,
        }
        self.data = {
            "interface_id": 2,
            "interface": {
                "h_bonds": {
                    "n_bonds": 2,
                    "bonds": [{**self.bond, "dist": 3.98765}, self.bond],
                },
                "salt_bridges": {
                    "n_bonds": 1,
                    "bonds": [{**self.bond, "dist": 4.32167}],
                },
                "ss_bonds": {
                    "n_bonds": 1,
                    "bonds": [
                        {
                            **self.bond,
                            "auth_comp_id_1": "CYS",
                            "auth_comp_id_2": "CYS",
                            "auth_atom_id_1": "SG",
                            "auth_atom_id_2": "SG",
                            "dist": 2.04123,
                        }
                    ],
                },
                "cov_bonds": {
                    "n_bonds": 1,
                    "bonds": [{**self.bond, "dist": 1.54321}],
                },
                "other_bonds": {"bonds": [{"unused": True}]},
            },
        }
        self.expected_tables = {
            "h_bonds": [{**self.expected_bond, "dist": 3.988}, self.expected_bond],
            "salt_bridges": [{**self.expected_bond, "dist": 4.322}],
            "cov_bonds": [{**self.expected_bond, "dist": 1.543}],
            "ss_bonds": [
                {
                    **self.expected_bond,
                    "label_comp_id_1": "CYS",
                    "label_comp_id_2": "CYS",
                    "auth_atom_id_1": "SG",
                    "auth_atom_id_2": "SG",
                    "dist": 2.041,
                }
            ],
        }

    def _write_interface(self, directory, data, compressed=False):
        os.makedirs(directory, exist_ok=True)
        extension = ".json.gz" if compressed else ".json"
        opener = gzip.open if compressed else open
        path = os.path.join(directory, f"interface_{data['interface_id']}" + extension)
        with opener(path, "wt") as f:
            json.dump(data, f)

    def test_parse_bond_tables(self):
        """Write four lists per interface with correct fields, order, and rounding."""
        for compressed in (False, True):
            with (
                self.subTest(compressed=compressed),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                inputs = os.path.join(tmpdir, "interfaces")
                output = os.path.join(tmpdir, "tables", "bonds")
                self._write_interface(inputs, self.data, compressed)
                self._write_interface(
                    inputs, {"interface_id": 10, "interface": {}}, compressed
                )
                extension = ".json.gz" if compressed else ".json"
                for name in ("summary", "interface_notes", "interface_2_h_bonds"):
                    with open(os.path.join(inputs, name + extension), "w") as f:
                        f.write("not JSON")
                os.mkdir(os.path.join(inputs, "interface_3" + extension))
                other_extension = ".json" if compressed else ".json.gz"
                with open(
                    os.path.join(inputs, "interface_99" + other_extension), "w"
                ) as f:
                    f.write("not JSON")

                PostProcessInterfacingBondsTable(inputs, output, compressed).parse()

                self.assertEqual(
                    set(os.listdir(output)),
                    {
                        f"interface_{interface_id}_{bond_type}" + extension
                        for interface_id in (2, 10)
                        for bond_type in self.expected_tables
                    },
                )
                opener = gzip.open if compressed else open
                for interface_id in (2, 10):
                    for bond_type, expected in self.expected_tables.items():
                        with opener(
                            os.path.join(
                                output,
                                f"interface_{interface_id}_{bond_type}" + extension,
                            ),
                            "rt",
                        ) as f:
                            self.assertEqual(
                                json.load(f), expected if interface_id == 2 else []
                            )

    def test_empty_bond_sections(self):
        """Missing or null bond lists still produce all four empty table files."""
        for section in (None, {}, {"n_bonds": 0}, {"bonds": None}, {"bonds": []}):
            with self.subTest(section=section), tempfile.TemporaryDirectory() as tmpdir:
                inputs = os.path.join(tmpdir, "interfaces")
                output = os.path.join(tmpdir, "bonds")
                data = {
                    "interface_id": 1,
                    "interface": {
                        bond_type: section for bond_type in self.expected_tables
                    },
                }
                self._write_interface(inputs, data)

                PostProcessInterfacingBondsTable(inputs, output).parse()

                self.assertEqual(len(os.listdir(output)), 4)
                for bond_type in self.expected_tables:
                    with open(
                        os.path.join(output, f"interface_1_{bond_type}.json")
                    ) as f:
                        self.assertEqual(json.load(f), [])

    def test_missing_or_empty_directory(self):
        """Match the other directory parsers by skipping absent or empty inputs."""
        for missing in (False, True):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as tmpdir:
                inputs = os.path.join(tmpdir, "interfaces")
                output = os.path.join(tmpdir, "bonds")
                if not missing:
                    os.makedirs(inputs)

                with self.assertLogs(
                    "pisa_utils.post_process_parsers", level="WARNING"
                ):
                    PostProcessInterfacingBondsTable(inputs, output).parse()

                self.assertFalse(os.path.exists(output))

    def test_invalid_bonds(self):
        """Validate all four tables before writing any file for the interface."""
        for bond_type in self.expected_tables:
            with (
                self.subTest(bond_type=bond_type),
                tempfile.TemporaryDirectory() as tmpdir,
            ):
                data = deepcopy(self.data)
                data["interface"][bond_type]["bonds"][0]["dist"] = "invalid"
                inputs = os.path.join(tmpdir, "interfaces")
                output = os.path.join(tmpdir, "bonds")
                self._write_interface(inputs, data)

                with self.assertRaises(ValidationError):
                    PostProcessInterfacingBondsTable(inputs, output).parse()

                self.assertEqual(os.listdir(output), [])

    def test_mismatched_interface_id(self):
        """Reject inconsistent IDs without writing files for the wrong interface."""
        with tempfile.TemporaryDirectory() as tmpdir:
            inputs = os.path.join(tmpdir, "interfaces")
            output = os.path.join(tmpdir, "bonds")
            self._write_interface(inputs, self.data)
            os.rename(
                os.path.join(inputs, "interface_2.json"),
                os.path.join(inputs, "interface_3.json"),
            )

            with self.assertRaisesRegex(ValueError, "Interface ID mismatch"):
                PostProcessInterfacingBondsTable(inputs, output).parse()

            self.assertEqual(os.listdir(output), [])

    def test_parse_real_interfaces(self):
        """Project real PISA bond data and keep empty disulfide/covalent tables."""
        inputs = "tests/data/expected_output/interfaces/3hax_interfaces_multi"
        with tempfile.TemporaryDirectory() as output:
            PostProcessInterfacingBondsTable(inputs, output).parse()

            self.assertEqual(len(os.listdir(output)), 12)
            with open(os.path.join(output, "interface_1_h_bonds.json")) as f:
                bonds = json.load(f)
            self.assertEqual(len(bonds), 35)
            self.assertEqual(
                bonds[0],
                {
                    "auth_asym_id_1": "E",
                    "label_comp_id_1": "G",
                    "auth_seq_id_1": 41,
                    "auth_atom_id_1": "N2",
                    "auth_asym_id_2": "A",
                    "label_comp_id_2": "GLU",
                    "auth_seq_id_2": 64,
                    "auth_atom_id_2": "OE2",
                    "dist": 2.844,
                },
            )
            with open(os.path.join(output, "interface_2_salt_bridges.json")) as f:
                self.assertEqual(len(json.load(f)), 18)
            for interface_id in (1, 2, 3):
                for bond_type in ("ss_bonds", "cov_bonds"):
                    with open(
                        os.path.join(
                            output, f"interface_{interface_id}_{bond_type}.json"
                        )
                    ) as f:
                        self.assertEqual(json.load(f), [])
