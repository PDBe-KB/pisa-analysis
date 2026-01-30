from unittest import TestCase

from pisa_utils.models.data_models import Complex


class TestComplexModel(TestCase):
    maxDiff = None

    def test_no_asm_set_no_asu_complex(self):
        """
        Test parsing a complex model where no assembly results are available due to
        overlapping assemblies. This is a known edge case of PISA output for some
        entries.
        """

        data = {
            "name": "TEST",
            "status": "Ok",
            "status_description": (
                "Results not available due to the following: \n Overlapping "
            ),
            "status_note": (
                "If this message is issued in error, report it \n to Eugene"
            ),
        }

        complex_model = Complex(**data)

        self.assertDictEqual(
            complex_model.model_dump(),
            {
                "session_name": "TEST",
                "status": "Ok",
                "status_description": (
                    "Results not available due to the following: Overlapping"
                ),
                "status_note": (
                    "If this message is issued in error, report it to Eugene"
                ),
                "all_chains_at_identity": None,
                "assessment": None,
                "asu_complex": None,
                "multimeric_state": None,
                "n_interfaces": None,
                "n_pqs_sets": None,
                "pqs_sets": [],
                "resolution": None,
            },
        )
