from unittest import TestCase

from pisa_utils.models.data_models import InterfaceSummary, Residue


class TestResidueModel(TestCase):
    maxDiff = None

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


class TestInterfaceSummary(TestCase):
    maxDiff = None

    def test_minimal_non_parsing(self):
        data = {
            "n_interfaces": 1,
            "interface_types": [
                {
                    "int_type": 1,
                    "interfaces": [
                        {
                            "interface_id": 1,
                            "auth_asym_id_1": "A",
                            "int_natoms_1": 100,
                            "int_nres_1": 1000,
                            "auth_asym_id_2": "B",
                            "int_natoms_2": 200,
                            "int_nres_2": 2000,
                            "int_area": 150000.0,
                            "int_solv_energy": -5000.0,
                            "pvalue": 0.05,
                            "css": 0.123,
                            "complex_keys_with_interface": [1, 2, 3],
                        },
                    ],
                },
            ],
        }

        interface_summary = InterfaceSummary(**data)
        self.assertDictEqual(
            interface_summary.model_dump(),
            {
                "n_interfaces": 1,
                "interface_types": [
                    {
                        "int_type": 1,
                        "interfaces": [
                            {
                                "interface_id": 1,
                                "auth_asym_id_1": "A",
                                "int_natoms_1": 100,
                                "int_nres_1": 1000,
                                "auth_asym_id_2": "B",
                                "int_natoms_2": 200,
                                "int_nres_2": 2000,
                                "int_area": 150000.0,
                                "int_solv_energy": -5000.0,
                                "pvalue": 0.05,
                                "css": 0.123,
                                "complex_keys_with_interface": [1, 2, 3],
                            },
                        ],
                    },
                ],
            },
        )

    def test_interface_with_string_numbers(self):
        data = {
            "n_interfaces": "2",
            "interface_types": [
                {
                    "int_type": "1",
                    "interfaces": [
                        {
                            "interface_id": "2",
                            "auth_asym_id_1": "C",
                            "int_natoms_1": "150",
                            "int_nres_1": "1500",
                            "auth_asym_id_2": "D",
                            "int_natoms_2": "250",
                            "int_nres_2": "2500",
                            "int_area": "175000.5678",
                            "int_solv_energy": "-6000.9876",
                            "pvalue": "0.01",
                            "css": "0.4567",
                            "complex_keys_with_interface": ["4", "5", "6"],
                        },
                    ],
                },
            ],
        }

        interface_summary = InterfaceSummary(**data)
        self.assertDictEqual(
            interface_summary.model_dump(),
            {
                "n_interfaces": 2,
                "interface_types": [
                    {
                        "int_type": 1,
                        "interfaces": [
                            {
                                "interface_id": 2,
                                "auth_asym_id_1": "C",
                                "int_natoms_1": 150,
                                "int_nres_1": 1500,
                                "auth_asym_id_2": "D",
                                "int_natoms_2": 250,
                                "int_nres_2": 2500,
                                "int_area": 175000.6,
                                "int_solv_energy": -6001.0,
                                "pvalue": 0.01,
                                "css": 0.457,
                                "complex_keys_with_interface": [4, 5, 6],
                            },
                        ],
                    },
                ],
            },
        )

    def test_missing_optional_fields(self):
        data = {
            "n_interfaces": 1,
            "interface_types": [
                {
                    "int_type": 1,
                    "interfaces": [
                        {
                            "interface_id": 1,
                            "auth_asym_id_1": "A",
                            "int_natoms_1": 100,
                            "int_nres_1": 1000,
                            "auth_asym_id_2": "B",
                            "int_natoms_2": 200,
                            "int_nres_2": 2000,
                            "int_area": 150000.0,
                            "int_solv_energy": -5000.0,
                            "pvalue": 0.05,
                            # missing optional fields
                        },
                    ],
                },
            ],
        }

        interface_summary = InterfaceSummary(**data)
        self.assertDictEqual(
            interface_summary.model_dump(),
            {
                "n_interfaces": 1,
                "interface_types": [
                    {
                        "int_type": 1,
                        "interfaces": [
                            {
                                "interface_id": 1,
                                "auth_asym_id_1": "A",
                                "int_natoms_1": 100,
                                "int_nres_1": 1000,
                                "auth_asym_id_2": "B",
                                "int_natoms_2": 200,
                                "int_nres_2": 2000,
                                "int_area": 150000.0,
                                "int_solv_energy": -5000.0,
                                "pvalue": 0.05,
                                "css": None,
                                "complex_keys_with_interface": None,
                            },
                        ],
                    },
                ],
            },
        )

    def test_multiple_interface_types(self):
        data = {
            "n_interfaces": 2,
            "interface_types": [
                {
                    "int_type": 1,
                    "interfaces": [
                        {
                            "interface_id": 1,
                            "auth_asym_id_1": "A",
                            "int_natoms_1": 100,
                            "int_nres_1": 1000,
                            "auth_asym_id_2": "B",
                            "int_natoms_2": 200,
                            "int_nres_2": 2000,
                            "int_area": 150000.0,
                            "int_solv_energy": -5000.0,
                            "pvalue": 0.05,
                            "css": 0.123,
                            "complex_keys_with_interface": [1, 2, 3],
                        },
                        {
                            "interface_id": 3,
                            "auth_asym_id_1": "A",
                            "int_natoms_1": 120,
                            "int_nres_1": 1100,
                            "auth_asym_id_2": "C",
                            "int_natoms_2": 220,
                            "int_nres_2": 2100,
                            "int_area": 160000.0,
                            "int_solv_energy": -5500.0,
                            "pvalue": 0.03,
                            "css": 0.234,
                            "complex_keys_with_interface": [7, 8, 9],
                        },
                    ],
                },
                {
                    "int_type": 2,
                    "interfaces": [
                        {
                            "interface_id": 2,
                            "auth_asym_id_1": "C",
                            "int_natoms_1": 150,
                            "int_nres_1": 1500,
                            "auth_asym_id_2": "D",
                            "int_natoms_2": 250,
                            "int_nres_2": 2500,
                            "int_area": 175000.0,
                            "int_solv_energy": -6000.0,
                            "pvalue": 0.01,
                            "css": 0.456,
                            "complex_keys_with_interface": [4, 5, 6],
                        },
                        {
                            "interface_id": 4,
                            "auth_asym_id_1": "D",
                            "int_natoms_1": 180,
                            "int_nres_1": 1800,
                            "auth_asym_id_2": "E",
                            "int_natoms_2": 280,
                            "int_nres_2": 2800,
                            "int_area": 185000.0,
                            "int_solv_energy": -6500.0,
                            "pvalue": 0.02,
                            "css": 0.567,
                            "complex_keys_with_interface": [10, 11, 12],
                        },
                    ],
                },
            ],
        }

        interface_summary = InterfaceSummary(**data)
        self.assertDictEqual(
            interface_summary.model_dump(),
            data,
        )
