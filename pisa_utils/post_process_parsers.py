from abc import ABC, abstractmethod
import json
import logging
import os
from typing import Optional

from pisa_utils.models.post_process_models import (
    InterfaceDetails,
    InterfaceDetailsComponent,
    InterfaceDetailsList,
    PQSSetRow,
    ComplexTableRow,
    ComplexTable,
    PISAAnalysisType,
)
from pisa_utils.file_io import save_json, open_compressed

LOGGER = logging.getLogger(__name__)


class PostProcessor(ABC):
    def __init__(
        self, input_json_path: str, output_json_path: str, compressed: bool = False
    ):
        self.input_json_path = input_json_path
        self.output_json_path = output_json_path
        self.compressed = compressed

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
        with open_compressed(self.input_json_path, compressed=self.compressed) as f:
            self.data = json.load(f)


class PostProcessComplexTable(PostProcessor):
    def parse(self):
        super().parse()

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


class PostProcessInterfaceDetailsList(PostProcessor):
    def __init__(
        self,
        path_interfaces: str,
        path_monomers_json: str,
        output_json_path: str,
        compressed: bool = False,
    ):
        super().__init__(
            input_json_path=path_monomers_json,
            output_json_path=output_json_path,
            compressed=compressed,
        )
        self.path_interfaces = path_interfaces

    def _extract_monomer_data(
        self, component_id: str, interface_auth_asym_id: Optional[str] = None
    ) -> dict:
        """
        Extract monomer (component) data for a given component ID

        :param component_id: ID of the component to extract data for. Not auth_asym_id
        :type component_id: str
        :param interface_auth_asym_id: auth_asym_id from the interface.json file. File
            provenance differs from monomers_extended.json, so used for sanity checks.
            Defaults to None
        :type interface_auth_asym_id: str, optional
        :raises ValueError: If no monomer data is found
        :return: Component (monomer) data for the specified component ID
        :rtype: dict
        """

        monomer_data_filtered = next(
            (
                item
                for item in self.data["components"]
                if item["component_id"] == component_id
            ),
            None,
        )

        # Sanity checks
        if not monomer_data_filtered:
            LOGGER.error(f"No monomer data found for component ID: {component_id}")
            raise ValueError(f"No monomer data found for component ID: {component_id}")

        if (
            interface_auth_asym_id
            and monomer_data_filtered["auth_asym_id"] != interface_auth_asym_id
        ):
            raise ValueError(
                f"auth_asym_id mismatch: interface says {interface_auth_asym_id}, "
                f"but monomer data is {monomer_data_filtered['auth_asym_id']}"
            )

        return monomer_data_filtered

    def parse(self):
        super().parse()
        LOGGER.info("Generating interface details file for all interfaces")

        output = []

        if not os.path.isdir(self.path_interfaces):
            LOGGER.warning(
                f"Interface JSON directory not found: {self.path_interfaces}. "
                "No interface details JSON will be created."
            )
            return

        ext = ".json.gz" if self.compressed else ".json"
        interface_files = [
            f
            for f in os.listdir(self.path_interfaces)
            if f.startswith("interface_") and f.endswith(ext)
            and os.path.isfile(os.path.join(self.path_interfaces, f))
        ]

        if not interface_files:
            LOGGER.warning(
                "No interface JSON files found in directory: "
                f"{self.path_interfaces}. No interface details JSON will be created."
            )
            return

        # Iterate over interfaces JSONs
        for interface_json_file in interface_files:
            LOGGER.debug(f"Processing interface JSON file: {interface_json_file}")

            with open_compressed(
                os.path.join(self.path_interfaces, interface_json_file),
                compressed=self.compressed,
            ) as f:
                interface_data = json.load(f)

            # Interface level data
            interface_id = interface_data["interface_id"]
            interface = interface_data["interface"]
            int_type = interface["int_type"]
            css = interface["css"]

            components: list[InterfaceDetailsComponent] = []

            for molecule in interface.get("molecules", []):
                auth_asym_id = molecule["auth_asym_id"]
                component_id = molecule["component_id"]

                LOGGER.debug(f"Processing molecule with component_id: {component_id}")

                monomer_data_filtered = self._extract_monomer_data(
                    component_id=component_id, interface_auth_asym_id=auth_asym_id
                )

                # Component-level data
                component = InterfaceDetailsComponent(
                    auth_asym_id=auth_asym_id,
                    molecule_class=molecule["molecule_class"],
                    symmetry_operation=molecule["symmetry_operation"],
                    symmetry_id=molecule["symmetry_id"],
                    int_natoms=molecule["int_natoms"],
                    surface_natoms=monomer_data_filtered["n_surface_atoms"],
                    total_atoms=monomer_data_filtered["total_atoms"],
                    int_nres=molecule["int_nres"],
                    surface_nres=monomer_data_filtered["n_surface_residues"],
                    total_residues=monomer_data_filtered["total_residues"],
                    int_area=molecule["int_area"],
                    surface_area=monomer_data_filtered["asa"],
                    solv_energy=monomer_data_filtered.get("solv_energy", None),
                    int_solv_energy=molecule["int_solv_energy"],
                    pvalue=molecule["pvalue"],
                )
                components.append(component)

            # Add interface details to ordered list
            interface_details = InterfaceDetails(
                interface_id=interface_id,
                int_type=int_type,
                css=css,
                components=components,
            )
            output.append(interface_details)

        # Order output list by interface ID for consistency
        output.sort(key=lambda x: x.interface_id)

        interface_details_list = InterfaceDetailsList(output).model_dump()

        self.save(self.output_json_path, interface_details_list)
