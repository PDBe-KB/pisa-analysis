from abc import ABC, abstractmethod
import json
import logging
import os
import re
from typing import Optional

from pisa_utils.models.post_process_models import (
    InterfaceDetails,
    InterfaceDetailsComponent,
    ComplexTableRow,
    ComplexTable,
    PISAAnalysisType,
    InterfaceTable,
    InterfaceComponentResidueTable,
    InterfacingBondRow,
    InterfacingBondsTable,
    ResidueRow,
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
    """Write a flat table of ASU and PQS complexes with per-row provenance."""

    def parse(self):
        super().parse()

        output = []

        # Add ASU complex
        asu_complex = self.data.get("asu_complex", {}).get("complex")
        if asu_complex:
            LOGGER.info("Processing ASU complex")
            asu_complex_row = ComplexTableRow(
                pqs_set_id=None,
                pisa_analysis_type=PISAAnalysisType.ASU,
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
            output.append(asu_complex_row)

        # Add PQS sets
        for pqs_set in self.data.get("pqs_sets", []):
            pqs_set_id = pqs_set["pqs_set_id"]
            complexes = pqs_set.get("complexes", [])
            LOGGER.info(
                f"Processing PQS set ID: {pqs_set_id} with {len(complexes)} complexes"
            )
            output.extend(
                ComplexTableRow(
                    pqs_set_id=int(pqs_set_id),
                    pisa_analysis_type=PISAAnalysisType.PQS,
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
            )

        complex_table = ComplexTable(output).model_dump()

        self.save(self.output_json_path, complex_table)


class PostProcessInterfaceDetailsList(PostProcessor):
    """Join interface and monomer data and write one InterfaceDetails file per ID.

    Output files are named interface_<id>.json, or .json.gz when compressed.
    """

    def __init__(
        self,
        path_interfaces: str,
        path_monomers_json: str,
        output_json_dir: str,
        compressed: bool = False,
    ):
        super().__init__(
            input_json_path=path_monomers_json,
            output_json_path=output_json_dir,
            compressed=compressed,
        )
        self.path_interfaces = path_interfaces
        self.output_json_dir = output_json_dir

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
        LOGGER.info("Generating individual interface details files")

        if not os.path.isdir(self.path_interfaces):
            LOGGER.warning(
                f"Interface JSON directory not found: {self.path_interfaces}. "
                "No interface details JSON will be created."
            )
            return

        ext = ".json.gz" if self.compressed else ".json"
        interface_files = [
            f
            for f in sorted(os.listdir(self.path_interfaces))
            if re.fullmatch(r"interface_\d+" + re.escape(ext), f)
            and os.path.isfile(os.path.join(self.path_interfaces, f))
        ]

        if not interface_files:
            LOGGER.warning(
                "No interface JSON files found in directory: "
                f"{self.path_interfaces}. No interface details JSON will be created."
            )
            return

        os.makedirs(self.output_json_dir, exist_ok=True)

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

            interface_details = InterfaceDetails(
                interface_id=interface_id,
                int_type=int_type,
                css=css,
                components=components,
            )
            self.save(
                os.path.join(
                    self.output_json_dir,
                    f"interface_{interface_details.interface_id}" + ext,
                ),
                interface_details.model_dump(),
            )


class PostProcessInterfaceTable(PostProcessor):
    data: InterfaceTable

    def __init__(
        self,
        input_summary_interfaces_path: str,
        path_interface_jsons: str,
        output_json_path: str,
        compressed: bool = False,
    ):
        super().__init__(input_summary_interfaces_path, output_json_path, compressed)

        self.path_interface_jsons = path_interface_jsons

    def parse(self):
        """
        Read the grouped interface_summary.json file and write a flat table.
        Bond counts are read from individual interface files in path_interface_jsons.
        Each row includes int_type and is ordered by int_type, then interface_id.
        """
        super().parse()

        LOGGER.debug("Parsing interface table")

        # Fetch counts of specific bond types.
        rows = []
        extension = ".json.gz" if self.compressed else ".json"
        for interface_type in self.data["interface_types"]:
            for interface in interface_type["interfaces"]:
                row = {**interface, "int_type": interface_type["int_type"]}
                interface_id = interface["interface_id"]
                path_interface_json = os.path.join(
                    self.path_interface_jsons, f"interface_{interface_id}" + extension
                )

                with open_compressed(
                    path_interface_json, compressed=self.compressed
                ) as f:
                    interface_obj = json.load(f)
                    interface_data = interface_obj.get("interface", {})

                    row["n_h_bonds"] = interface_data.get("h_bonds", {}).get(
                        "n_bonds", 0
                    )
                    row["n_salt_bridges"] = interface_data.get("salt_bridges", {}).get(
                        "n_bonds", 0
                    )
                    row["n_ss_bonds"] = interface_data.get("ss_bonds", {}).get(
                        "n_bonds", 0
                    )
                rows.append(row)

        self.data = InterfaceTable.model_validate(rows)

        # Sort interfaces by type and ID.
        self.data.root.sort(key=lambda row: (row.int_type, row.interface_id))

        self.save(self.output_json_path, self.data.model_dump())


class PostProcessInterfaceComponentResidueTable(PostProcessor):
    """Write one component residue table for each interface JSON in a directory.

    Input files must be named interface_<id>.json (or .json.gz when compressed).
    Output files are named interface_<id>_residue_table.json, with the same
    compression setting. Molecule and residue order are preserved.
    """

    def __init__(
        self,
        path_interfaces: str,
        output_json_dir: str,
        compressed: bool = False,
    ):
        super().__init__(path_interfaces, output_json_dir, compressed)
        self.path_interfaces = path_interfaces
        self.output_json_dir = output_json_dir

    def _parse_residues(self, molecule: dict) -> list[ResidueRow]:
        """Project molecule residues onto the table's fields and derive metrics."""
        return [
            ResidueRow(
                auth_asym_id=molecule["auth_asym_id"],
                label_comp_id=residue["auth_comp_id"],
                auth_seq_id=residue["auth_seq_id"],
                bonds=residue.get("bonds"),
                asa=residue["asa"],
                bsa=residue["bsa"],
                solv_energy=residue["solv_energy"],
                # Required inputs; the model derives these from rounded ASA/BSA.
                is_interfacing_residue=False,
                is_solvent_accessible_residue=False,
                is_inaccessible_residue=False,
            )
            for residue in molecule["residues"]["residues"]
        ]

    def parse(self) -> None:
        if not os.path.isdir(self.path_interfaces):
            LOGGER.warning(
                f"Interface JSON directory not found: {self.path_interfaces}. "
                "No residue table JSON files will be created."
            )
            return

        extension = ".json.gz" if self.compressed else ".json"
        interface_files = [
            filename
            for filename in sorted(os.listdir(self.path_interfaces))
            if re.fullmatch(r"interface_\d+" + re.escape(extension), filename)
            and os.path.isfile(os.path.join(self.path_interfaces, filename))
        ]
        if not interface_files:
            LOGGER.warning(
                f"No interface JSON files found in directory: {self.path_interfaces}. "
                "No residue table JSON files will be created."
            )
            return

        os.makedirs(self.output_json_dir, exist_ok=True)

        for filename in interface_files:
            with open_compressed(
                os.path.join(self.path_interfaces, filename),
                compressed=self.compressed,
            ) as f:
                interface_data = json.load(f)

            file_interface_id = int(filename[len("interface_") : -len(extension)])
            interface_id = int(interface_data["interface_id"])
            if interface_id != file_interface_id:
                raise ValueError(
                    f"Interface ID mismatch: filename has {file_interface_id}, JSON has {interface_id}"
                )
            molecules = interface_data["interface"]["molecules"]
            if len(molecules) != 2:
                raise ValueError(
                    f"Interface {interface_id} must contain exactly two molecules; "
                    f"found {len(molecules)}"
                )

            table = InterfaceComponentResidueTable(
                structure_1_residues=self._parse_residues(molecules[0]),
                structure_2_residues=self._parse_residues(molecules[1]),
            )
            self.save(
                os.path.join(
                    self.output_json_dir,
                    f"interface_{interface_id}_residue_table" + extension,
                ),
                table.model_dump(),
            )


class PostProcessInterfacingBondsTable(PostProcessor):
    """Write hydrogen bond, salt bridge, disulfide, and covalent tables per interface.

    Read interface_<id>.json files and write interface_<id>_<bond_type>.json
    for h_bonds, salt_bridges, ss_bonds, and cov_bonds. Empty types produce empty lists.
    With compression enabled, both input and output filenames end in .json.gz.
    """

    BOND_TYPES = ("h_bonds", "salt_bridges", "ss_bonds", "cov_bonds")

    def __init__(
        self,
        path_interfaces: str,
        output_json_dir: str,
        compressed: bool = False,
    ):
        super().__init__(path_interfaces, output_json_dir, compressed)
        self.path_interfaces = path_interfaces
        self.output_json_dir = output_json_dir

    def _parse_bonds(self, bonds: list[dict]) -> InterfacingBondsTable:
        """Keep the table fields and let the model round bond distances."""
        return InterfacingBondsTable(
            [
                InterfacingBondRow(
                    auth_asym_id_1=bond["auth_asym_id_1"],
                    label_comp_id_1=bond["auth_comp_id_1"],
                    auth_seq_id_1=bond["auth_seq_id_1"],
                    auth_atom_id_1=bond["auth_atom_id_1"],
                    dist=bond["dist"],
                    auth_asym_id_2=bond["auth_asym_id_2"],
                    label_comp_id_2=bond["auth_comp_id_2"],
                    auth_seq_id_2=bond["auth_seq_id_2"],
                    auth_atom_id_2=bond["auth_atom_id_2"],
                )
                for bond in bonds
            ]
        )

    def parse(self) -> None:
        if not os.path.isdir(self.path_interfaces):
            LOGGER.warning(
                f"Interface JSON directory not found: {self.path_interfaces}. "
                "No bond table JSON files will be created."
            )
            return

        extension = ".json.gz" if self.compressed else ".json"
        interface_files = [
            filename
            for filename in sorted(os.listdir(self.path_interfaces))
            if re.fullmatch(r"interface_\d+" + re.escape(extension), filename)
            and os.path.isfile(os.path.join(self.path_interfaces, filename))
        ]
        if not interface_files:
            LOGGER.warning(
                f"No interface JSON files found in directory: {self.path_interfaces}. "
                "No bond table JSON files will be created."
            )
            return

        os.makedirs(self.output_json_dir, exist_ok=True)

        for filename in interface_files:
            with open_compressed(
                os.path.join(self.path_interfaces, filename),
                compressed=self.compressed,
            ) as f:
                interface_data = json.load(f)

            file_interface_id = int(filename[len("interface_") : -len(extension)])
            interface_id = int(interface_data["interface_id"])
            if interface_id != file_interface_id:
                raise ValueError(
                    f"Interface ID mismatch: filename has {file_interface_id}, "
                    f"JSON has {interface_id}"
                )

            interface = interface_data["interface"]
            tables = {}
            for bond_type in self.BOND_TYPES:
                bond_data = interface.get(bond_type) or {}
                tables[bond_type] = self._parse_bonds(bond_data.get("bonds") or [])

            for bond_type, table in tables.items():
                self.save(
                    os.path.join(
                        self.output_json_dir,
                        f"interface_{interface_id}_{bond_type}" + extension,
                    ),
                    table.model_dump(),
                )
