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
    LOGGER.warning(f"Iterating over {len(error_details)} pydantic validation errors...")

    for i, error in enumerate(error_details):
        # auth_asym_id (chain ID) not assigned by user
        if error["input"] is None and error["loc"][-1] in CHAIN_ID_SYNONYMS:
            LOGGER.warning(
                f"Pydantic error {i}: Missing chain ID detected. Details: {error}"
            )
            raise MissingChainIDError()

        # Future, bespoke validation error handling...
