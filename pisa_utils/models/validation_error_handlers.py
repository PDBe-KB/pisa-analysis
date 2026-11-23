from enum import Enum

from pydantic import ValidationError
from pydantic_core import ErrorDetails
from logging import getLogger

LOGGER = getLogger(__name__)

CHAIN_ID_SYNONYMS = {"chain-1", "chain-2", "chain_id", "monomer_1", "monomer_2"}


class MissingChainIDError(Exception):
    """Custom exception for missing chain ID in the input data."""

    def __init__(self, message: str = None):

        if not message:
            message = (
                "File contains unassigned chain IDs. Please make sure all polymer and "
                "ligand entities in the file have assigned alphanumeric chain IDs."
            )

        super().__init__(message)


def trigger_helpful_validation_error(error_details: list[ErrorDetails]):
    """
    Custom validation error handler to check for specific validation errors and raise
    bespoke exceptions for API to handle accordingly.

    :param error_details: List of validation error details from pydantic validation.
    :type error_details: list[ErrorDetails]
    :raises MissingChainIDError: If any validation error indicates missing chain ID
        assignment in the input data.
    """

    for i, error in enumerate(error_details):
        # Error due to auth_asym_id (chain ID) not assigned by user
        if error["input"] is None and error["loc"][-1] in CHAIN_ID_SYNONYMS:
            LOGGER.warning(
                f"Pydantic error {i}: Missing chain ID detected. Details: {error}"
            )
            raise MissingChainIDError()

        LOGGER.warning(f"Pydantic error {i}: Validation error details: {error}")
        raise ValidationError.from_exception_data(
            "Unhandled pydantic validation errors", error_details
        )


class ErrorLabels(str, Enum):
    # Raised by failures from PISA binary
    PISA_NO_RESULTS = "no_results"
    PISA_STOPPED = "pisa_stopped"
    PISA_CONFIGURATION_ERROR = "configuration_error"
    PISA_SESSION_DIR_FAILURE = "session_dir_failure"
    PISA_NO_SRS = "no_srs"
    PISA_NO_MOL_REF_INDEX = "no_mol_ref_index"
    PISA_NO_MOL_REF_DATA = "no_mol_ref_data"
    PISA_NO_AGENTS_DATA = "no_agents_data"
    PISA_NOT_ANALYSED = "not_analysed_or_file_not_found"
    PISA_NO_ATOM_DATA_FOUND = "no_atom_data_found"
    PISA_NO_COORDINATE_DATA = "no_coordinate_data"
    PISA_NO_MACROMOLECULES = "no_macromolecules"
    PISA_NOT_READY_FOR_CONTACTS = "not_ready_for_contacts"
    PISA_NO_SPACE_GROUP = "no_space_group"
    PISA_NO_RCSB_SPACE_GROUP = "no_rcsb_space_group"
    PISA_NO_SYM_OPS = "no_sym_ops"
    PISA_NO_SYM_OP = "no_sym_op"
    PISA_SURFACE_CALCS_FAILED = "surface_calcs_failed"
    PISA_BAD_COORDINATE_FILE = "bad_coordinate_file"
    PISA_CANT_WRITE_RESULTS = "cant_write_results"
    PISA_CANT_READ_RESULTS = "cant_read_results"
    PISA_SRS_FAILED = "srs_failed"
    PISA_INTERFACE_CALC_FAILURE = "interface_calc_failure"
    PISA_NO_RCSB_DATA = "no_rcsb_data"
    PISA_NO_RCSB_ASSIGNMENT = "no_rcsb_assignment"
    PISA_NO_ASSEMBLY_PARAMETERS = "no_assembly_parameters"
    PISA_NO_INTERFACE_STATS = "no_interface_stats"
    PISA_UNKNOWN_INTERFACE_TYPE = "unknown_interface_type"
    PISA_WRONG_REF_PARAMETER = "wrong_ref_parameter"
    PISA_INCOMPLETE_STAT_SET = "incomplete_stat_set"
    PISA_UNKNOWN_INTF_STAT_ERROR = "unknown_intf_stat_error"
    PISA_SESSION_DOESNT_EXIST = "session_doesnt_exist"
    PISA_NO_SESSION_RESULTS = "no_session_results"
    PISA_CANT_WRITE_FILE = "cant_write_file"
    PISA_SESSION_ERASED = "session_erased"
    PISA_MONOMER_NO_OUT_OF_RANGE = "monomer_no_out_of_range"
    PISA_INTERFACE_NO_OUT_OF_RANGE = "interface_no_out_of_range"
    PISA_ASSEMBLY_NO_OUT_OF_RANGE = "assembly_no_out_of_range"
    PISA_CANT_WRITE_RASMOL_DATA = "cant_write_rasmol_data"
    PISA_CANT_WRITE_MG_DATA = "cant_write_mg_data"
    PISA_NO_DIMER_FOUND = "no_dimer_found"
    PISA_UNKNOWN = "unknown_pisa_error"

    # From parsers defined by pisa-analysis
    PARSING = "parsing_error"

    # Fallback
    UNKNOWN = "unknown_error"


KNOWN_PISA_ERROR_CODES = {
    # Relevant error codes taken from:
    # https://github.com/PDBe-KB/pisa/blob/6d6f7d09ea5bf0f3f8a1602dd56fa896e6aca401/pisalib/pisa_types.h#L47
    -2: ErrorLabels.PISA_STOPPED,
    -1: ErrorLabels.PISA_NO_RESULTS,
    3: ErrorLabels.PISA_CONFIGURATION_ERROR,
    4: ErrorLabels.PISA_SESSION_DIR_FAILURE,
    5: ErrorLabels.PISA_NO_SRS,
    6: ErrorLabels.PISA_NO_MOL_REF_INDEX,
    7: ErrorLabels.PISA_NO_MOL_REF_DATA,
    8: ErrorLabels.PISA_NO_AGENTS_DATA,
    9: ErrorLabels.PISA_NOT_ANALYSED,
    10: ErrorLabels.PISA_NO_ATOM_DATA_FOUND,
    11: ErrorLabels.PISA_NO_COORDINATE_DATA,
    12: ErrorLabels.PISA_NO_MACROMOLECULES,
    13: ErrorLabels.PISA_NOT_READY_FOR_CONTACTS,
    14: ErrorLabels.PISA_NO_SPACE_GROUP,
    15: ErrorLabels.PISA_NO_RCSB_SPACE_GROUP,
    16: ErrorLabels.PISA_NO_SYM_OPS,
    17: ErrorLabels.PISA_NO_SYM_OP,
    18: ErrorLabels.PISA_SURFACE_CALCS_FAILED,
    19: ErrorLabels.PISA_BAD_COORDINATE_FILE,
    20: ErrorLabels.PISA_CANT_WRITE_RESULTS,
    21: ErrorLabels.PISA_CANT_READ_RESULTS,
    22: ErrorLabels.PISA_SRS_FAILED,
    23: ErrorLabels.PISA_INTERFACE_CALC_FAILURE,
    24: ErrorLabels.PISA_NO_RCSB_DATA,
    25: ErrorLabels.PISA_NO_RCSB_ASSIGNMENT,
    26: ErrorLabels.PISA_NO_ASSEMBLY_PARAMETERS,
    27: ErrorLabels.PISA_NO_INTERFACE_STATS,
    28: ErrorLabels.PISA_UNKNOWN_INTERFACE_TYPE,
    29: ErrorLabels.PISA_WRONG_REF_PARAMETER,
    30: ErrorLabels.PISA_INCOMPLETE_STAT_SET,
    31: ErrorLabels.PISA_UNKNOWN_INTF_STAT_ERROR,
    32: ErrorLabels.PISA_SESSION_DOESNT_EXIST,
    33: ErrorLabels.PISA_NO_SESSION_RESULTS,
    34: ErrorLabels.PISA_CANT_WRITE_FILE,
    35: ErrorLabels.PISA_SESSION_ERASED,
    36: ErrorLabels.PISA_MONOMER_NO_OUT_OF_RANGE,
    37: ErrorLabels.PISA_INTERFACE_NO_OUT_OF_RANGE,
    38: ErrorLabels.PISA_ASSEMBLY_NO_OUT_OF_RANGE,
    39: ErrorLabels.PISA_CANT_WRITE_RASMOL_DATA,
    40: ErrorLabels.PISA_CANT_WRITE_MG_DATA,
    41: ErrorLabels.PISA_NO_DIMER_FOUND,
    # Parsing error code for pisa-analysis parsers
    999: ErrorLabels.PARSING,
}
