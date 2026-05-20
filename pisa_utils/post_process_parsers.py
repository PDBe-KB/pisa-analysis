from abc import ABC, abstractmethod
import json
import logging

from pisa_utils.models.post_process_models import (
    PQSSetRow,
    ComplexTableRow,
    ComplexTable,
    PISAAnalysisType,
)
from pisa_utils.file_io import save_json

LOGGER = logging.getLogger(__name__)


class PostProcessor(ABC):
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

        save_json(data, path_json, compressed=self.compressed)

        LOGGER.info(f"JSON file written successfully: {path_json}")

    @abstractmethod
    def parse(self):
        pass


class PostProcessComplexTable(PostProcessor):
    def parse(self):
        output = []

        # Add ASU complex
        asu_complex = self.data.get("asu_complex", {}).get("complex")
        if asu_complex:
            LOGGER.info("Processing ASU complex")
            asu_complex_row = ComplexTableRow(
                complex_key=asu_complex["complex_key"],
                formula=asu_complex["formula"],
                composition=asu_complex["composition"],
                asa=asu_complex["asa"],
                bsa=asu_complex["bsa"],
                int_energy=asu_complex["int_energy"],
                diss_energy=asu_complex["diss_energy"],
                entropy=asu_complex["entropy"],
                mmsize=asu_complex["mmsize"],
                n_uc=asu_complex["n_uc"],
                symmetry_number=asu_complex["symmetry_number"],
                n_interfaces=asu_complex.get("interfaces", {}).get("n_interfaces"),
            )
            asu_row = PQSSetRow(
                pisa_analysis_type=PISAAnalysisType.ASU, complexes=[asu_complex_row]
            )
            output.append(asu_row)

        # Add PQS sets
        for pqs_set in self.data.get("pqs_sets", []):
            pqs_set_id = pqs_set["pqs_set_id"]
            complexes = pqs_set.get("complexes", [])
            LOGGER.info(
                f"Processing PQS set ID: {pqs_set_id} with {len(complexes)} complexes"
            )
            pqs_set_row = PQSSetRow(
                pqs_set_id=int(pqs_set_id),
                pisa_analysis_type=PISAAnalysisType.PQS,
                complexes=[
                    ComplexTableRow(
                        complex_key=complex["complex_key"],
                        formula=complex["formula"],
                        composition=complex["composition"],
                        asa=complex["asa"],
                        bsa=complex["bsa"],
                        int_energy=complex["int_energy"],
                        diss_energy=complex["diss_energy"],
                        entropy=complex["entropy"],
                        mmsize=complex["mmsize"],
                        n_uc=complex["n_uc"],
                        symmetry_number=complex["symmetry_number"],
                        n_interfaces=complex.get("interfaces", {}).get("n_interfaces"),
                    )
                    for complex in complexes
                ],
            )
            output.append(pqs_set_row)

        complex_table = ComplexTable(output).model_dump()

        self.save(self.output_json_path, complex_table)
