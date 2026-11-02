from unittest import TestCase
from copy import deepcopy
import os
import json

from pydantic import ValidationError

from pisa_utils.models.post_process_models import (
    ComplexTable,
    InterfaceDetailsList,
    InterfaceTable,
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

        ComplexTable.model_validate(data)

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
        self.assertIsNone(complex_table[0].complexes[0].formula)
        self.assertEqual(complex_table[0].complexes[0].n_interfaces, 0)

        # Check PQS set values
        self.assertEqual(complex_table[1].pisa_analysis_type.value, "PQS set")
        self.assertEqual(complex_table[1].pqs_set_id, 1)
        self.assertIsNone(complex_table[1].complexes[0].formula)
        self.assertEqual(complex_table[1].complexes[0].n_interfaces, 0)


class TestInterfaceDetailsList(PostProcessTestCase):
    def test_interface_details_list(self):
        """
        Simple model validation. Optional values omitted from the input.
        """

        input_json_path = os.path.join(
            self.path_input_stubs, "interface_details_list.json"
        )
        data = self._open_json(input_json_path)

        interface_details_list = InterfaceDetailsList.model_validate(data).root

        self.assertIsNone(interface_details_list[0].components[0].solv_energy)
        self.assertEqual(interface_details_list[0].components[0].pvalue, 0.973)
        self.assertEqual(
            interface_details_list[0].components[0].int_solv_energy, -20.833
        )
        self.assertEqual(interface_details_list[0].components[0].surface_area, 11475.9)
        self.assertEqual(interface_details_list[0].css, 0.661)


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
        """
        interface_table = InterfaceTable.model_validate(self.data)

        expected_data = deepcopy(self.data)
        expected_data["interface_types"][0]["interfaces"][0].update(
            int_area=1617.4,
            int_solv_energy=-23.1,
            pvalue=0.661,
            css=0.382,
        )
        self.assertDictEqual(interface_table.model_dump(), expected_data)

    def test_interface_table_aliases_and_optional_fields(self):
        """
        Accept input aliases with optional fields omitted or explicitly null.

        Resolve chain_id and int_solv_en to the model's field names, keep CSS
        and complex_keys_with_interface as None, and omit aliases from output.
        """
        for optional_fields in ("omitted", "null"):
            with self.subTest(optional_fields=optional_fields):
                data = deepcopy(self.data)
                row = data["interface_types"][0]["interfaces"][0]
                row["chain_id"] = row.pop("auth_asym_id_2")
                row["int_solv_en"] = row.pop("int_solv_energy")
                for field in ("css", "complex_keys_with_interface"):
                    if optional_fields == "omitted":
                        row.pop(field)
                    else:
                        row[field] = None

                interface_table = InterfaceTable.model_validate(data)
                interface = interface_table.interface_types[0].interfaces[0]

                self.assertEqual(interface.auth_asym_id_2, "A")
                self.assertEqual(interface.int_solv_energy, -23.1)
                self.assertIsNone(interface.css)
                self.assertIsNone(interface.complex_keys_with_interface)
                self.assertNotIn("chain_id", interface.model_dump())
                self.assertNotIn("int_solv_en", interface.model_dump())

    def test_interface_table_rejects_invalid_data(self):
        """
        Reject missing required fields and unexpected fields at each level.

        Check table, interface type group, and interface row validation errors,
        asserting that each reports the expected field location and error type.
        """
        cases = (
            ((), "n_interfaces", "missing"),
            (("interface_types", 0), "int_type", "missing"),
            (("interface_types", 0, "interfaces", 0), "interface_id", "missing"),
            ((), "unexpected", "extra_forbidden"),
            (("interface_types", 0), "unexpected", "extra_forbidden"),
            (("interface_types", 0, "interfaces", 0), "unexpected", "extra_forbidden"),
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
