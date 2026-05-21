from unittest import TestCase
import os
import json
from pisa_utils.models.post_process_models import ComplexTable


class TestComplexTable(TestCase):
    def setUp(self):
        super().setUp()
        self.base_path = "tests/data/"
        self.path_input_stubs = os.path.join(
            self.base_path, "mock_data/model_inputs/post_processor_models/"
        )

    def _open_json(self, file_path):
        with open(file_path, "r") as f:
            return json.load(f)

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
