from unittest import TestCase

from pisa_utils.models.data_models import Residue


class TestResidueModel(TestCase):
    def test_with_aliased_fields(self):
        data = {
            "ser_no": 1,
            "name": "ALA",
            "seq_num": 1,
            "label_seq_num": 1,
            "ins_code": "A",
            "bonds": "C",
            "asa": 100.0,
            "bsa": 100.0,
            "solv_en": -10.0,
        }

        residue = Residue(**data)
        self.assertDictEqual(
            residue.model_dump(),
            {
                "residue_serial_number": 1,
                "auth_comp_id": "ALA",
                "auth_seq_id": 1,
                "label_seq_id": 1,
                "ins_code": "A",
                "bonds": "C",
                "asa": 100.0,
                "bsa": 100.0,
                "solv_energy": -10.0,
            },
        )

    def test_with_attribute_fields(self):
        data = {
            "residue_serial_number": 2,
            "auth_comp_id": "GLY",
            "auth_seq_id": 2,
            "label_seq_id": 2,
            "ins_code": "B",
            "bonds": "N",
            "asa": 150.0,
            "bsa": 150.0,
            "solv_energy": -15.0,
        }

        residue = Residue(**data)
        self.assertDictEqual(
            residue.model_dump(),
            {
                "residue_serial_number": 2,
                "auth_comp_id": "GLY",
                "auth_seq_id": 2,
                "label_seq_id": 2,
                "ins_code": "B",
                "bonds": "N",
                "asa": 150.0,
                "bsa": 150.0,
                "solv_energy": -15.0,
            },
        )

    def test_missing_optional_fields(self):
        data = {
            "residue_serial_number": 1,
            "auth_comp_id": "ALA",
            "auth_seq_id": 1,
            "asa": 100.0,
            "bsa": 100.0,
            "solv_energy": 100.0,
        }

        residue = Residue(**data)
        self.assertDictEqual(
            residue.model_dump(),
            {
                "residue_serial_number": 1,
                "auth_comp_id": "ALA",
                "auth_seq_id": 1,
                "label_seq_id": None,
                "ins_code": None,
                "bonds": None,
                "asa": 100.0,
                "bsa": 100.0,
                "solv_energy": 100.0,
            },
        )

    def test_string_to_number_conversion(self):
        data = {
            "ser_no": "3",
            "name": "SER",
            "seq_num": "3",
            "label_seq_num": "3",
            "ins_code": "C",
            "bonds": "O",
            "asa": "200.001234",
            "bsa": "200.001234",
            "solv_en": "-20.001234",
        }

        residue = Residue(**data)
        self.assertDictEqual(
            residue.model_dump(),
            {
                "residue_serial_number": 3,
                "auth_comp_id": "SER",
                "auth_seq_id": 3,
                "label_seq_id": 3,
                "ins_code": "C",
                "bonds": "O",
                "asa": 200.0,
                "bsa": 200.0,
                "solv_energy": -20.001,
            },
        )
