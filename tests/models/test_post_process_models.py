from unittest import TestCase

from pisa_utils.models.post_process_models import ComplexTable


class TestComplexTable(TestCase):
    def test_complex_table_asu_and_multi_pqs_sets(self):
        data = [
            {
                "pqs_set_id": None,
                "pisa_analysis_type": "Asymmetric unit",
                "complexes": [
                    {
                        "complex_key": 0,
                        "formula": "ABCDEa(2)b(2)de(5)f",
                        "composition": "ACDEF[PGE](2)[EDO](2)[PG4][MG](5)[FHU]",
                        "asa": 30301.0,
                        "bsa": 14879.2,
                        "int_energy": -78.5,
                        "diss_energy": 27.891,
                        "entropy": 11.7,
                        "mmsize": 5,
                        "n_uc": 0,
                        "symmetry_number": 1,
                        "n_interfaces": 23,
                    }
                ],
            },
            {
                "pqs_set_id": 1,
                "pisa_analysis_type": "PQS set",
                "complexes": [
                    {
                        "complex_key": 1,
                        "formula": "ABCDEa(2)b(2)cde(5)f",
                        "composition": "ACDEF[PGE](2)[EDO](2)[ZN][PG4][MG](5)[FHU]",
                        "asa": 30770.9,
                        "bsa": 14507.2,
                        "int_energy": -94.4,
                        "diss_energy": 27.836,
                        "entropy": 11.8,
                        "mmsize": 5,
                        "n_uc": 4,
                        "symmetry_number": 1,
                        "n_interfaces": 25,
                    }
                ],
            },
            {
                "pqs_set_id": 2,
                "pisa_analysis_type": "PQS set",
                "complexes": [
                    {
                        "complex_key": 2,
                        "formula": "BCDa(2)bde(5)",
                        "composition": "CDE[PGE](2)[EDO][PG4][MG](5)",
                        "asa": 20209.7,
                        "bsa": 4495.9,
                        "int_energy": -43.5,
                        "diss_energy": 4.838,
                        "entropy": 11.3,
                        "mmsize": 3,
                        "n_uc": 4,
                        "symmetry_number": 1,
                        "n_interfaces": 15,
                    },
                    {
                        "complex_key": 3,
                        "formula": "AEbcf",
                        "composition": "AF[EDO][ZN][FHU]",
                        "asa": 18173.4,
                        "bsa": 2399.0,
                        "int_energy": -16.4,
                        "diss_energy": 16.211,
                        "entropy": 10.6,
                        "mmsize": 2,
                        "n_uc": 4,
                        "symmetry_number": 1,
                        "n_interfaces": 5,
                    },
                ],
            },
            {
                "pqs_set_id": 3,
                "pisa_analysis_type": "PQS set",
                "complexes": [
                    {
                        "complex_key": 4,
                        "formula": "CDEabde(5)f",
                        "composition": "DEF[PGE][EDO][PG4][MG](5)[FHU]",
                        "asa": 19082.6,
                        "bsa": 3962.4,
                        "int_energy": -33.9,
                        "diss_energy": 1.398,
                        "entropy": 10.6,
                        "mmsize": 3,
                        "n_uc": 4,
                        "symmetry_number": 1,
                        "n_interfaces": 13,
                    },
                    {
                        "complex_key": 5,
                        "formula": "ABabc",
                        "composition": "AC[PGE][EDO][ZN]",
                        "asa": 19200.6,
                        "bsa": 3032.5,
                        "int_energy": -11.7,
                        "diss_energy": 12.095,
                        "entropy": 11.4,
                        "mmsize": 2,
                        "n_uc": 4,
                        "symmetry_number": 1,
                        "n_interfaces": 4,
                    },
                    {
                        "complex_key": 6,
                        "formula": "ABEabcf",
                        "composition": "ACF[PGE][EDO][ZN][FHU]",
                        "asa": 20774.3,
                        "bsa": 5245.9,
                        "int_energy": -27.3,
                        "diss_energy": 12.056,
                        "entropy": 11.4,
                        "mmsize": 3,
                        "n_uc": 4,
                        "symmetry_number": 1,
                        "n_interfaces": 7,
                    },
                ],
            },
        ]

        ComplexTable.model_validate(data)

    def test_complex_table_no_pqs_sets(self):
        data = [
            {
                "pqs_set_id": None,
                "pisa_analysis_type": "Asymmetric unit",
                "complexes": [
                    {
                        "complex_key": 0,
                        "formula": "ABCDEa(2)b(2)de(5)f",
                        "composition": "ACDEF[PGE](2)[EDO](2)[PG4][MG](5)[FHU]",
                        "asa": 30301.0,
                        "bsa": 14879.2,
                        "int_energy": -78.5,
                        "diss_energy": 27.891,
                        "entropy": 11.7,
                        "mmsize": 5,
                        "n_uc": 0,
                        "symmetry_number": 1,
                        "n_interfaces": 23,
                    }
                ],
            }
        ]

        ComplexTable.model_validate(data)

    def test_complex_table_only_mandatory_fields(self):
        data = [
            {
                "pqs_set_id": None,
                "pisa_analysis_type": "Asymmetric unit",
                "complexes": [
                    {
                        "complex_key": 0,
                        "composition": "ACDEF[PGE](2)[EDO](2)[PG4][MG](5)[FHU]",
                        "asa": 30301.0,
                        "bsa": 14879.2,
                        "int_energy": -78.5,
                        "diss_energy": 27.891,
                        "entropy": 11.7,
                        "mmsize": 5,
                        "n_uc": 0,
                        "symmetry_number": 1,
                    }
                ],
            },
            {
                "pqs_set_id": 1,
                "complexes": [
                    {
                        "complex_key": 1,
                        "composition": "ACDEF[PGE](2)[EDO](2)[ZN][PG4][MG](5)[FHU]",
                        "asa": 30770.9,
                        "bsa": 14507.2,
                        "int_energy": -94.4,
                        "diss_energy": 27.836,
                        "entropy": 11.8,
                        "mmsize": 5,
                        "n_uc": 4,
                        "symmetry_number": 1,
                    }
                ],
            },
        ]

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
