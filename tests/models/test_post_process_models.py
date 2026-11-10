from unittest import TestCase
from copy import deepcopy
import os
import json

from pydantic import ValidationError

from pisa_utils.models.post_process_models import (
    ComplexTable,
    InterfaceDetails,
    InterfaceTable,
    ResidueRow,
)


class PostProcessTestCase(TestCase):
    def setUp(self):
        super().setUp()
        self.base_path = "tests/data/"
        self.path_input_stubs = os.path.join(
            self.base_path, "mock_data/model_inputs/post_processor_models/"
        )

    def _open_json(self, file_path):
        with open(file_path, "r") as f:
            return json.load(f)


class TestComplexTable(PostProcessTestCase):
    def test_complex_table_asu_and_multi_pqs_sets(self):
        """
        Test that the ComplexTable model can handle both ASU and multiple PQS sets,
        with all optional fields provided
        """

        input_json_path = os.path.join(self.path_input_stubs, "multi_pqs.json")
        data = self._open_json(input_json_path)

        complex_table = ComplexTable.model_validate(data)

        self.assertEqual(complex_table.model_dump(mode="json"), data)
        self.assertEqual(
            [row.pqs_set_id for row in complex_table.root],
            [None, 1, 2, 2, 3, 3, 3],
        )

    def test_complex_table_no_pqs_sets(self):
        """
        Test for only asymmetric unit (ASU) data, no PQS sets found situation.
        """

        input_json_path = os.path.join(self.path_input_stubs, "only_asu.json")
        data = self._open_json(input_json_path)

        ComplexTable.model_validate(data)

    def test_complex_table_only_mandatory_fields(self):
        input_json_path = os.path.join(
            self.path_input_stubs, "optional_fields_omitted.json"
        )
        data = self._open_json(input_json_path)

        complex_table = ComplexTable.model_validate(data).root

        # Check ASU complex values
        self.assertEqual(complex_table[0].pisa_analysis_type.value, "Asymmetric unit")
        self.assertIsNone(complex_table[0].pqs_set_id)
        self.assertIsNone(complex_table[0].formula)
        self.assertEqual(complex_table[0].n_interfaces, 0)

        # Check PQS set values
        self.assertEqual(complex_table[1].pisa_analysis_type.value, "PQS set")
        self.assertEqual(complex_table[1].pqs_set_id, 1)
        self.assertIsNone(complex_table[1].formula)
        self.assertEqual(complex_table[1].n_interfaces, 0)


class TestInterfaceDetails(PostProcessTestCase):
    def test_interface_details(self):
        """
        Simple model validation. Optional values omitted from the input.
        """

        details = []
        for interface_id in (1, 2, 3):
            input_json_path = os.path.join(
                self.path_input_stubs,
                "interface_details",
                f"interface_{interface_id}.json",
            )
            data = self._open_json(input_json_path)
            detail = InterfaceDetails.model_validate(data)
            self.assertEqual(detail.interface_id, interface_id)
            details.append(detail)

        self.assertIsNone(details[0].components[0].solv_energy)
        self.assertEqual(details[0].components[0].pvalue, 0.973)
        self.assertEqual(details[0].components[0].int_solv_energy, -20.833)
        self.assertEqual(details[0].components[0].surface_area, 11475.9)
        self.assertEqual(details[0].css, 0.661)


class TestInterfaceTable(PostProcessTestCase):
    def setUp(self):
        super().setUp()
        input_json_path = os.path.join(self.path_input_stubs, "interface_table.json")
        self.data = self._open_json(input_json_path)

    def test_interface_table_rounding(self):
        """
        Validate a complete table and check its serialized contents.

        Round interface area and solvation energy to one decimal place, and
        pvalue and CSS to three, preserving all other fields and the structure.
        Keep explicit hydrogen bond, salt bridge, and disulfide bond counts.
        """
        interface_table = InterfaceTable.model_validate(self.data)

        expected_data = deepcopy(self.data)
        expected_data[0].update(
            int_area=1617.4,
            int_solv_energy=-23.1,
            pvalue=0.661,
            css=0.382,
        )
        self.assertListEqual(interface_table.model_dump(), expected_data)

    def test_interface_table_optional_fields(self):
        """
        Accept canonical field names with optional fields omitted or null.

        Default missing bond counts to zero and other optional fields to None;
        preserve explicit null values for all of them.
        """
        for optional_fields in ("omitted", "null"):
            with self.subTest(optional_fields=optional_fields):
                data = deepcopy(self.data)
                row = data[0]
                bond_fields = ("n_h_bonds", "n_salt_bridges", "n_ss_bonds")
                for field in ("css", "complex_keys_with_interface", *bond_fields):
                    if optional_fields == "omitted":
                        row.pop(field)
                    else:
                        row[field] = None

                interface_table = InterfaceTable.model_validate(data)
                interface = interface_table.root[0]

                self.assertEqual(interface.auth_asym_id_2, "A")
                self.assertEqual(interface.int_solv_energy, -23.1)
                self.assertIsNone(interface.css)
                self.assertIsNone(interface.complex_keys_with_interface)
                expected_count = 0 if optional_fields == "omitted" else None
                for field in bond_fields:
                    self.assertEqual(getattr(interface, field), expected_count)
                    self.assertEqual(interface.model_dump()[field], expected_count)
                self.assertNotIn("chain_id", interface.model_dump())
                self.assertNotIn("int_solv_en", interface.model_dump())

    def test_interface_table_rejects_removed_aliases(self):
        """Reject old aliases alone and alongside their canonical field names."""
        for field, alias in (
            ("auth_asym_id_2", "chain_id"),
            ("int_solv_energy", "int_solv_en"),
        ):
            for keep_canonical in (False, True):
                with self.subTest(alias=alias, keep_canonical=keep_canonical):
                    data = deepcopy(self.data)
                    data[0][alias] = data[0][field]
                    if not keep_canonical:
                        data[0].pop(field)

                    with self.assertRaises(ValidationError) as context:
                        InterfaceTable.model_validate(data)

                    expected_errors = {((0, alias), "extra_forbidden")}
                    if not keep_canonical:
                        expected_errors.add(((0, field), "missing"))
                    self.assertEqual(
                        {
                            (error["loc"], error["type"])
                            for error in context.exception.errors()
                        },
                        expected_errors,
                    )

    def test_interface_table_rejects_invalid_data(self):
        """
        Reject missing required fields and unexpected fields on each row.

        Check that each error identifies the row index, field, and error type.
        """
        cases = (
            ((0,), "int_type", "missing"),
            ((0,), "interface_id", "missing"),
            ((0,), "unexpected", "extra_forbidden"),
        )
        for path, field, error_type in cases:
            with self.subTest(path=path, field=field, error_type=error_type):
                data = deepcopy(self.data)
                target = data
                for key in path:
                    target = target[key]
                if error_type == "missing":
                    target.pop(field)
                else:
                    target[field] = True

                with self.assertRaises(ValidationError) as context:
                    InterfaceTable.model_validate(data)

                errors = context.exception.errors()
                self.assertEqual(len(errors), 1)
                self.assertEqual(errors[0]["loc"], path + (field,))
                self.assertEqual(errors[0]["type"], error_type)

    def test_interface_table_rejects_invalid_bond_counts(self):
        """Reject nonnumeric bond counts and identify the invalid row field."""
        for field in ("n_h_bonds", "n_salt_bridges", "n_ss_bonds"):
            with self.subTest(field=field):
                data = deepcopy(self.data)
                data[0][field] = "invalid"

                with self.assertRaises(ValidationError) as context:
                    InterfaceTable.model_validate(data)

                errors = context.exception.errors()
                self.assertEqual(len(errors), 1)
                self.assertEqual(errors[0]["loc"], (0, field))
                self.assertEqual(errors[0]["type"], "int_parsing")

    def test_interface_table_rejects_grouped_input(self):
        """The table's root must be a list rather than the old grouped object."""
        with self.assertRaises(ValidationError) as context:
            InterfaceTable.model_validate({"n_interfaces": 0, "interface_types": []})

        errors = context.exception.errors()
        self.assertEqual(len(errors), 1)
        self.assertEqual(errors[0]["loc"], ())
        self.assertEqual(errors[0]["type"], "list_type")


class TestResidueRow(PostProcessTestCase):
    def setUp(self):
        super().setUp()
        self.data = {
            "auth_asym_id": "A",
            "label_comp_id": "ARG",
            "auth_seq_id": 11,
            "asa": 158.712,
            "bsa": 123.567,
            "solv_energy": -1.23456,
            "is_interfacing_residue": False,
            "is_solvent_accessible_residue": False,
            "is_inaccessible_residue": False,
        }

    def test_residue_row_uses_solv_energy(self):
        """Validate and serialize the canonical energy field with rounding."""
        residue = ResidueRow.model_validate(self.data)

        self.assertEqual(residue.solv_energy, -1.235)
        self.assertEqual(residue.model_dump()["solv_energy"], -1.235)
        self.assertNotIn("solv_en", residue.model_dump())

    def test_residue_row_rejects_removed_alias(self):
        """Reject solv_en alone and alongside solv_energy."""
        for keep_canonical in (False, True):
            with self.subTest(keep_canonical=keep_canonical):
                data = deepcopy(self.data)
                data["solv_en"] = data["solv_energy"]
                if not keep_canonical:
                    data.pop("solv_energy")

                with self.assertRaises(ValidationError) as context:
                    ResidueRow.model_validate(data)

                expected_errors = {(("solv_en",), "extra_forbidden")}
                if not keep_canonical:
                    expected_errors.add((("solv_energy",), "missing"))
                self.assertEqual(
                    {
                        (error["loc"], error["type"])
                        for error in context.exception.errors()
                    },
                    expected_errors,
                )
