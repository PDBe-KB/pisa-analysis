from unittest import TestCase
from pydantic_core import ErrorDetails

from pisa_utils.models.validation_error_handlers import (
    ErrorLabels,
    KNOWN_PISA_ERROR_CODES,
    MissingChainIDError,
    trigger_helpful_validation_error,
)


def test_known_pisa_error_codes():
    expected = {
        -2: "pisa_stopped",
        -1: "no_results",
        3: "configuration_error",
        4: "session_dir_failure",
        5: "no_srs",
        6: "no_mol_ref_index",
        7: "no_mol_ref_data",
        8: "no_agents_data",
        9: "not_analysed_or_file_not_found",
        10: "no_atom_data_found",
        11: "no_coordinate_data",
        12: "no_macromolecules",
        13: "not_ready_for_contacts",
        14: "no_space_group",
        15: "no_rcsb_space_group",
        16: "no_sym_ops",
        17: "no_sym_op",
        18: "surface_calcs_failed",
        19: "bad_coordinate_file",
        20: "cant_write_results",
        21: "cant_read_results",
        22: "srs_failed",
        23: "interface_calc_failure",
        24: "no_rcsb_data",
        25: "no_rcsb_assignment",
        26: "no_assembly_parameters",
        27: "no_interface_stats",
        28: "unknown_interface_type",
        29: "wrong_ref_parameter",
        30: "incomplete_stat_set",
        31: "unknown_intf_stat_error",
        32: "session_doesnt_exist",
        33: "no_session_results",
        34: "cant_write_file",
        35: "session_erased",
        36: "monomer_no_out_of_range",
        37: "interface_no_out_of_range",
        38: "assembly_no_out_of_range",
        39: "cant_write_rasmol_data",
        40: "cant_write_mg_data",
        41: "no_dimer_found",
        999: "parsing_error",
    }

    assert all(
        isinstance(label, ErrorLabels) for label in KNOWN_PISA_ERROR_CODES.values()
    )
    assert {
        code: label.value for code, label in KNOWN_PISA_ERROR_CODES.items()
    } == expected


class TestTriggerHelpfulValidationError(TestCase):
    def test_trigger_missing_chain_in_interface_validation(self):
        error_details = [
            ErrorDetails(
                **{
                    "type": "string_type",
                    "loc": ("interface", "other-bonds", "bond", 0, "chain-2"),
                    "msg": "Input should be a valid string",
                    "input": None,
                    "url": "https://errors.pydantic.dev/2.12/v/string_type",
                }
            )
        ]

        with self.assertRaises(MissingChainIDError) as context:  # noqa: F841
            trigger_helpful_validation_error(error_details)
