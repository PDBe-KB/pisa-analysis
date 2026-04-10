from unittest import TestCase
from pydantic_core import ErrorDetails

from pisa_utils.models.validation_error_handlers import (
    MissingChainIDError,
    trigger_helpful_validation_error,
)


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
