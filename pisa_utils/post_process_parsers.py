from abc import ABC, abstractmethod
import json
import logging
import gzip

from pisa_utils.models.post_process_models import (
    PQSSetRow,
    ComplexTableRow,
    ComplexTable,
)

LOGGER = logging.getLogger(__name__)


class PostProcesser(ABC):
    def __init__(self, input_json_path, output_json_path, compressed=False):
        self.input_json_path = input_json_path
        self.output_json_path = output_json_path
        self.compressed = compressed

        with open(self.input_json_path, "r") as f:
            self.data = json.load(f)

    def save(self, path_json: str, data: dict) -> None:
        """
        Save data to JSON file, with optional gzip compression.

        :param path_json: Path to the output JSON file.
        :type path_json: str
        :param data: Data to save to JSON.
        :type data: dict
        """

        if self.compressed and not path_json.endswith(".gz"):
            path_json += ".gz"

        if self.compressed:
            with gzip.open(path_json, "wt") as json_file:
                json.dump(data, json_file, indent=4)
        else:
            with open(path_json, "w") as json_file:
                json.dump(data, json_file, indent=4)

        LOGGER.info(f"JSON file written successfully: {path_json}")

    @abstractmethod
    def parse(self):
        pass


class PostProcessComplexTable(PostProcesser):
    def parse(self):
        output = []

        for pqs_set in self.data.get("pqs_sets", []):
            pqs_set_id = pqs_set["pqs_set_id"]
            complexes = pqs_set.get("complexes", [])
            LOGGER.info(
                f"Processing PQS set ID: {pqs_set_id} with {len(complexes)} complexes"
            )
            pqs_set_row = PQSSetRow(
                pqs_set_id=int(pqs_set_id),
                complexes=[
                    ComplexTableRow(
                        complex_key=complex["complex_key"],
                        formula=complex["formula"],
                        composition=complex["composition"],
                        asa=complex["asa"],
                        bsa=complex["bsa"],
                        int_energy=complex["int_energy"],
                        diss_energy=complex["diss_energy"],
                    )
                    for complex in complexes
                ],
            )
            output.append(pqs_set_row)

        complex_table = ComplexTable(output).model_dump()

        self.save(self.output_json_path, complex_table)
