import json
import logging
import os
from abc import ABC, abstractmethod

import gemmi
import pandas as pd
import xmltodict

from pisa_utils.models.data_models import (
    Complex,
    ComplexExtended,
    Components,
    Interface,
    InterfaceExtended,
    InterfaceSummary,
)
from pisa_utils.utils import extract_ligand_contents, id_is_ligand, is_int_or_float

LOGGER = logging.getLogger(__name__)


def get_residue_mappings(structure: gemmi.Structure) -> dict:
    """
    Generate a mapping of auth_asym_id to residue ranges for polymers and ligands.

    :param structure: Gemmi Structure object.
    :type structure: gemmi.Structure
    :return: Mapping of auth_asym_id to residue ranges.
    :rtype: dict
    """
    residue_mappings = {}

    for chain in structure[0]:
        auth_asym_id = chain.name
        residue_mappings[auth_asym_id] = {
            "polymers": {},
            "ligands": {},
        }

        polymer = chain.get_polymer()
        res_seq_ids = [residue.seqid.num for residue in polymer]
        residue_mappings[auth_asym_id]["polymers"] = {
            "start": min(res_seq_ids),
            "end": max(res_seq_ids),
        }

        for ligand in chain.get_ligands():
            ccd_id = ligand.name
            auth_seq_id = ligand.seqid.num
            residue_mappings[auth_asym_id]["ligands"][ccd_id] = auth_seq_id

    return residue_mappings


class ConvertXMLToJSON(ABC):
    def __init__(
        self,
        path_xml: str,
        path_output: str,
        path_structure_file: str,
        data_type: str = None,
    ) -> None:
        """
        Base class for converting PISA-generated XML files to JSON files.

        :param path_xml: Path to the XML file.
        :type path_xml: str
        :param path_output: Path to the output JSON file or directory.
        :type path_output: str
        :param path_structure_file: Path to the input structure file (PDB or mmCIF).
        :type path_structure_file: str
        :param data_type: Possible options are defined by the inheriting object,
            defaults to None
        :type data_type: str, optional
        """
        self.path_xml = path_xml
        self.path_output = path_output
        self.path_structure_file = path_structure_file
        self.data_type = data_type

        # Instantiated during runtime
        self.df: pd.DataFrame

    def __str__(self):
        return (
            f"Convert{self.data_type}XMLToJSON: {self.path_xml} -> {self.path_output}"
        )

    @abstractmethod
    def parse(self) -> dict:
        """
        Convert any XML file to a JSON file.
        """

        with open(self.path_xml) as xml_file:
            # Parse XML to JSON
            LOGGER.info(f"Reading XML file: {self.path_xml}")
            pisa_data_xml = xmltodict.parse(xml_file.read())

        return pisa_data_xml

    def _file_is_pdb(self) -> bool:
        """
        Determine if the structure file is a PDB file.

        :return: True if PDB file, False otherwise.
        :rtype: bool
        """
        return self.path_structure_file.endswith(
            ".pdb"
        ) or self.path_structure_file.endswith(".ent")

    def _extract_mmcif_contents(self, path_mmcif: str) -> pd.DataFrame:
        """
        Extract relevant contents from mmCIF file into a pandas DataFrame.

        :param path_mmcif: Path to the mmCIF file.
        :type path_mmcif: str
        :return: DataFrame containing relevant mmCIF contents.
        :rtype: pd.DataFrame
        """
        block = gemmi.cif.read(path_mmcif).sole_block()

        headers = [
            "group_PDB",
            "label_asym_id",
            "label_seq_id",
            "auth_seq_id",
            "auth_comp_id",
            "auth_asym_id",
        ]
        table = block.find("_atom_site.", headers)

        return pd.DataFrame(
            table,
            columns=headers,
        ).drop_duplicates()

    def _extract_ligand_rows(
        self, ccd_id: str, auth_asym_id: str, auth_seq_id: int
    ) -> pd.DataFrame:
        """
        Extract rows corresponding to a specific ligand from the mmCIF DataFrame.

        :param ccd_id: CCD ID of the ligand.
        :type ccd_id: str
        :param auth_asym_id: auth_asym_id from the mmCIF file or chain ID in PDB file.
        :type auth_asym_id: str
        :param auth_seq_id: auth_seq_id from the mmCIF file.
        :type auth_seq_id: int
        :return: DataFrame containing rows for the specified ligand.
        :rtype: pd.DataFrame
        """

        return self.df.loc[
            (self.df["group_PDB"] == "HETATM")
            & (self.df["auth_comp_id"] == ccd_id)
            & (self.df["auth_asym_id"] == auth_asym_id)
            & (self.df["auth_seq_id"].astype(int) == auth_seq_id)
        ]

    def _extract_polymer_rows(self, auth_asym_id: str) -> pd.DataFrame:
        """
        Extract rows corresponding to a specific polymer chain from the mmCIF DataFrame.

        :param auth_asym_id: auth_asym_id from the mmCIF file or chain ID in PDB file.
        :type auth_asym_id: str
        :return: DataFrame containing rows for the specified polymer chain.
        :rtype: pd.DataFrame
        """
        return self.df.loc[
            (self.df["group_PDB"] == "ATOM") & (self.df["auth_asym_id"] == auth_asym_id)
        ]

    def _validate_ligand_label_seq_id(
        self,
        ligand_df: pd.DataFrame,
        auth_seq_id: str,
        auth_asym_id: str,
        ccd_id: str = None,
    ) -> str:
        """
        Extract and validate label_seq_id from ligand dataframe.

        :param ligand_df: DataFrame containing ligand rows.
        :type ligand_df: pd.DataFrame
        :param auth_seq_id: auth_seq_id as defined in the mmCIF file.
        :type auth_seq_id: str
        :param auth_asym_id: auth_asym_id from the mmCIF file or chain ID in PDB file.
        :type auth_asym_id: str
        :param ccd_id: CCD ID of the ligand, defaults to None
        :type ccd_id: str, optional
        :return: label_seq_id corresponding to the ligand.
        :rtype: str
        """

        label_seq_ids = ligand_df["label_seq_id"].unique()

        if len(label_seq_ids) == 0:
            LOGGER.error(
                f"Could not find label_seq_id for ligand {ccd_id} "
                f"in chain {auth_asym_id} at residue {auth_seq_id}"
            )
        elif len(label_seq_ids) > 1:
            LOGGER.error(
                f"Multiple label_seq_id found for ligand {ccd_id} "
                f"in chain {auth_asym_id} at residue {auth_seq_id}: "
                f"{label_seq_ids}"
            )
        else:
            label_seq_id = label_seq_ids[0]
            return label_seq_id if label_seq_id != "." else None

    def _validate_ligand_label_asym_id(
        self,
        ligand_df: pd.DataFrame,
        auth_seq_id: int,
        auth_asym_id: str,
        ccd_id: str = None,
    ) -> str:
        """
        Extract and validate label_asym_id from ligand dataframe.

        :param ligand_df: DataFrame containing ligand rows.
        :type ligand_df: pd.DataFrame
        :param auth_seq_id: auth_seq_id as defined in the mmCIF file.
        :type auth_seq_id: int
        :param auth_asym_id: auth_asym_id from the mmCIF file or chain ID in PDB file.
        :type auth_asym_id: str
        :param ccd_id: CCD ID of the ligand, defaults to None
        :type ccd_id: str, optional
        :return: label_asym_id corresponding to the ligand.
        :rtype: str
        """

        label_asym_id = ligand_df["label_asym_id"].unique()
        if len(label_asym_id) == 0:
            LOGGER.error(
                f"Could not find label_asym_id for ligand {ccd_id} "
                f"in chain {auth_asym_id} at residue {auth_seq_id}"
            )
        elif len(label_asym_id) > 1:
            LOGGER.error(
                f"Multiple label_asym_id found for ligand {ccd_id} "
                f"in chain {auth_asym_id} at residue {auth_seq_id}: "
                f"{label_asym_id}"
            )
        else:
            return label_asym_id[0]

    def _extract_ligand_ids(
        self,
        ligand_df: pd.DataFrame,
        auth_seq_id: int,
        auth_asym_id: str,
        ccd_id: str = None,
    ) -> tuple[str, int]:
        """
        Extract label_asym_id and label_seq_id from ligand dataframe.

        :param ligand_df: DataFrame containing ligand rows.
        :type ligand_df: pd.DataFrame
        :param auth_seq_id: auth_seq_id as defined in the mmCIF file.
        :type auth_seq_id: int
        :param auth_asym_id: auth_asym_id from the mmCIF file or chain ID in PDB file.
        :type auth_asym_id: str
        :param ccd_id: CCD ID of the ligand, defaults to None
        :type ccd_id: str, optional
        :return: label_asym_id and label_seq_id corresponding to the ligand.
        :rtype: tuple[str, int]
        """

        # Extract label_asym_id from dataframe
        label_asym_id = self._validate_ligand_label_asym_id(
            ligand_df, auth_seq_id, auth_asym_id, ccd_id
        )

        # Extract first and last label_seq_ids from dataframe
        label_seq_id = self._validate_ligand_label_seq_id(
            ligand_df, auth_seq_id, auth_asym_id, ccd_id
        )

        return label_asym_id, label_seq_id

    def _validate_polymer_ids(
        self, polymer_df: pd.DataFrame, auth_asym_id: str
    ) -> tuple[str, int, int, str, str]:
        """
        Extract and validate label_asym_id, auth_seq_id_start, auth_seq_id_end,
        label_seq_id_start, and label_seq_id_end from polymer dataframe.

        :param polymer_df: DataFrame containing polymer rows for one chain.
        :type polymer_df: pd.DataFrame
        :param auth_asym_id: auth_asym_id from the mmCIF file or chain ID in PDB file.
        :type auth_asym_id: str
        :return: label_asym_id, auth_seq_id_start, auth_seq_id_end, label_seq_id_start,
            label_seq_id_end.
        :rtype: tuple[str, int, int, str, str]
        """

        # Extract label_asym_id from dataframe
        label_asym_ids = polymer_df["label_asym_id"].unique()

        if len(label_asym_ids) == 0:
            LOGGER.error(
                f"Could not find label_asym_id for polymer " f"in chain {auth_asym_id}"
            )
        elif len(label_asym_ids) > 1:
            LOGGER.error(
                f"Multiple label_asym_id found for polymer "
                f"in chain {auth_asym_id}: {label_asym_ids}"
            )
        else:
            label_asym_id = label_asym_ids[0]

        # Extract first and last residue numbers from dataframe
        auth_seq_id_start = polymer_df["auth_seq_id"].astype(int).min()
        auth_seq_id_end = polymer_df["auth_seq_id"].astype(int).max()
        label_seq_id_start = polymer_df["label_seq_id"].astype(int).min()
        label_seq_id_end = polymer_df["label_seq_id"].astype(int).max()

        return (
            label_asym_id,
            auth_seq_id_start,
            auth_seq_id_end,
            label_seq_id_start,
            label_seq_id_end,
        )


class ConvertInterfaceXMLToJSONs(ConvertXMLToJSON):
    def __init__(
        self, path_xml: str, path_jsons: str, path_structure_file: str
    ) -> None:
        """
        Handles conversion of a PISA-generated interface XML into N separate JSONs,
        each containing the interface data for interface N_i.

        :param path_xml: Path to the interface XML file.
        :type path_xml: str
        :param path_jsons: Path to the output directory for interface JSON files.
        :type path_jsons: str
        """
        super().__init__(path_xml, path_jsons, path_structure_file, "Interface")

        self.bond_types = (
            "h-bonds",
            "salt-bridges",
            "ss-bonds",
            "cov-bonds",
            "other-bonds",
        )

    def parse(self) -> None:
        """
        Converts interface file to set of JSON files.
        """

        all_interface_data = super().parse()

        self.pdb_id = all_interface_data["pdb_entry"]["pdb_code"]
        self.status = all_interface_data["pdb_entry"]["status"]
        self.n_interfaces = all_interface_data["pdb_entry"]["n_interfaces"]

        if self.n_interfaces == "0":
            LOGGER.warning(
                f"No interfaces found in XML file: {self.path_xml}. "
                "No JSON files will be created."
            )
            return

        # Make interfaces dir
        os.makedirs(self.path_output, exist_ok=True)
        LOGGER.info(f"Created directory for interfaces: {self.path_output}")

        # Map auth_asym_id to label_asym_id from structure file
        if not self._file_is_pdb():
            self.df = self._extract_mmcif_contents(self.path_structure_file)
        else:
            self.structure = gemmi.read_structure(self.path_structure_file)

        # Extract interfaces to separate JSON files
        if self.n_interfaces == "1":
            self._parse_interface_from_dict(
                all_interface_data["pdb_entry"]["interface"]
            )
        else:
            for interface_data in all_interface_data["pdb_entry"]["interface"]:
                self._parse_interface_from_dict(interface_data)

    def _parse_interface_from_dict(self, interface: dict) -> None:
        """
        Converts a single interface dictionary to a JSON file.

        :param interface: Dictionary containing interface data from original XML.
        :type interface: dict
        """
        interface_id = interface["id"]
        del interface["id"]

        interface_output = {
            "interface_id": interface_id,
            "n_interfaces": self.n_interfaces,
            "status": self.status,
            "pdb_id": self.pdb_id,
            "interface": interface,
        }

        # Set label_asym_* to None if PDB file
        if self._file_is_pdb():
            residue_mappings = get_residue_mappings(self.structure)

            for bond_type in self.bond_types:
                bond_info = interface_output["interface"].get(bond_type, [])
                bonds = bond_info.get("bond", [])

                if not bonds:
                    continue

                for bond in bonds if isinstance(bonds, list) else [bonds]:
                    bond["label_asym_id-1"] = None
                    bond["label_asym_id-2"] = None
                    bond["label_seqnum-1"] = None
                    bond["label_seqnum-2"] = None

            for molecule in interface_output["interface"].get("molecule", []):
                if not molecule:
                    continue

                auth_asym_id = molecule["chain_id"]
                if auth_asym_id in residue_mappings:
                    molecule["auth_seq_id_start"] = residue_mappings[auth_asym_id][
                        "polymers"
                    ]["start"]
                    molecule["auth_seq_id_end"] = residue_mappings[auth_asym_id][
                        "polymers"
                    ]["end"]

                residues = molecule.get("residues", {}).get("residue", [])
                for residue in residues if isinstance(residues, list) else [residues]:
                    if not residue:
                        continue
                    residue["label_seq_num"] = None

        # Add missing label_asym_ids from mmCIF file
        else:
            molecules = interface_output["interface"].get("molecule", [])
            for molecule in molecules if isinstance(molecules, list) else [molecules]:
                auth_asym_id = None
                ccd_id = None
                auth_seq_id_start = None
                auth_seq_id_end = None

                # Handle formatting for ligands
                if molecule["class"] == "Ligand":
                    (
                        auth_asym_id,
                        ccd_id,
                        auth_seq_id_start,
                    ) = extract_ligand_contents(molecule["chain_id"])

                    auth_seq_id_end = auth_seq_id_start

                    ligand_df = self._extract_ligand_rows(
                        ccd_id=ccd_id,
                        auth_asym_id=auth_asym_id,
                        auth_seq_id=auth_seq_id_start,
                    )

                    # Extract label_asym_id from dataframe
                    label_asym_id, label_seq_id = self._extract_ligand_ids(
                        ligand_df, auth_seq_id_start, auth_asym_id, ccd_id
                    )

                    label_seq_id_start = label_seq_id
                    label_seq_id_end = label_seq_id

                # Handle formatting for polymers
                else:
                    auth_asym_id = molecule["chain_id"]

                    polymer_df = self._extract_polymer_rows(auth_asym_id)

                    (
                        label_asym_id,
                        auth_seq_id_start,
                        auth_seq_id_end,
                        label_seq_id_start,
                        label_seq_id_end,
                    ) = self._validate_polymer_ids(polymer_df, auth_asym_id)

                # Set values in molecule dict
                molecule["chain_id"] = auth_asym_id
                molecule["label_asym_id"] = label_asym_id
                molecule["ccd_id"] = ccd_id
                molecule["auth_seq_id_start"] = auth_seq_id_start
                molecule["auth_seq_id_end"] = auth_seq_id_end
                molecule["label_seq_id_start"] = label_seq_id_start
                molecule["label_seq_id_end"] = label_seq_id_end

        # Validate through pydantic model
        interface_output = Interface(**interface_output).model_dump()

        # Write interface JSON
        interface_json_path = os.path.join(
            self.path_output, f"interface_{interface_id}.json"
        )

        with open(interface_json_path, "w") as interface_json_file:
            json.dump(interface_output, interface_json_file, indent=4)
            LOGGER.info(f"Interface JSON written: {interface_json_path}")


class ConvertAssemblyXMLToJSON(ConvertXMLToJSON):
    def __init__(
        self,
        path_xml: str,
        path_json: str,
        path_structure_file: str,
        path_interface_jsons: str,
    ) -> None:
        """
        Converts a PISA-generated assembly XML into a single JSON file.

        :param path_xml: Path to the assembly XML file.
        :type path_xml: str
        :param path_json: Path to the output assembly JSON file.
        :type path_json: str
        """
        super().__init__(path_xml, path_json, path_structure_file, "Assembly")
        self.path_interface_jsons = path_interface_jsons

    def parse(self) -> None:
        """
        Convert assembly XML file to single JSON file.
        """

        assembly_data = super().parse()

        structure = gemmi.read_structure(self.path_structure_file)

        # Add resolution
        assembly_data["pisa_results"]["resolution"] = structure.resolution

        # Count the number of interfaces across all asm_sets
        total_interfaces = 0
        asm_sets = assembly_data["pisa_results"].get("asm_set", [])
        if not asm_sets:
            total_interfaces = (
                assembly_data["pisa_results"]
                .get("asu_complex", {})
                .get("assembly", {})
                .get("interfaces", {})
                .get("n_interfaces", 0)
            )
        else:
            if not isinstance(asm_sets, list):
                asm_sets = [asm_sets]
            interface_ids = set()

            for asm_set in asm_sets:
                assembly = asm_set.get("assembly")

                if not isinstance(assembly, list):
                    assembly = [assembly]

                for asm in assembly:
                    interfaces = asm.get("interfaces", {}).get("interface", [])

                    if not isinstance(interfaces, list):
                        interfaces = [interfaces]
                    for interface in interfaces:
                        interface_ids.add(interface["id"])

            total_interfaces = len(interface_ids)

        assembly_data["pisa_results"]["n_interfaces"] = total_interfaces

        # Map auth_asym_id to label_asym_id from mmCIF file
        if not self._file_is_pdb():
            self.df = self._extract_mmcif_contents(self.path_structure_file)

            if isinstance(assembly_data["pisa_results"].get("asm_set"), dict):
                assembly_data["pisa_results"]["asm_set"] = [
                    assembly_data["pisa_results"]["asm_set"]
                ]

            # Add missing label_*_ids from mmCIF file to pqs_sets
            for pqs_set in assembly_data["pisa_results"].get("asm_set", []):
                if not pqs_set:
                    continue

                if isinstance(pqs_set.get("assembly"), dict):
                    pqs_set["assembly"] = [pqs_set["assembly"]]

                for complex in pqs_set.get("assembly", []):
                    if not complex:
                        continue

                    if isinstance(complex.get("molecule"), dict):
                        complex["molecule"] = [complex["molecule"]]

                    # Add missing label_*_ids to molecules in complex
                    for molecule in complex.get("molecule", []):
                        if not molecule:
                            continue

                        # Ligand molecule
                        if id_is_ligand(molecule["chain_id"]):
                            (
                                auth_asym_id,
                                ccd_id,
                                auth_seq_id,
                            ) = extract_ligand_contents(molecule["chain_id"])

                            ligand_df = self._extract_ligand_rows(
                                ccd_id=ccd_id,
                                auth_asym_id=auth_asym_id,
                                auth_seq_id=auth_seq_id,
                            )

                            # Extract label_asym_id from dataframe
                            label_asym_id, label_seq_id = self._extract_ligand_ids(
                                ligand_df, auth_seq_id, auth_asym_id, ccd_id
                            )

                            molecule["label_asym_id"] = label_asym_id
                            molecule["label_seq_id_start"] = label_seq_id
                            molecule["label_seq_id_end"] = label_seq_id

                        # Polymer molecule
                        else:
                            auth_asym_id = molecule["chain_id"]
                            polymer_df = self._extract_polymer_rows(auth_asym_id)

                            # Extract label_asym_id from dataframe
                            (
                                label_asym_id,
                                auth_seq_id_start,
                                auth_seq_id_end,
                                label_seq_id_start,
                                label_seq_id_end,
                            ) = self._validate_polymer_ids(polymer_df, auth_asym_id)

                            molecule["label_asym_id"] = label_asym_id
                            molecule["auth_seq_id_start"] = auth_seq_id_start
                            molecule["auth_seq_id_end"] = auth_seq_id_end
                            molecule["label_seq_id_start"] = label_seq_id_start
                            molecule["label_seq_id_end"] = label_seq_id_end

            # Add missing label_*_ids to ASU complex
            asu_molecules = (
                assembly_data["pisa_results"]
                .get("asu_complex", {})
                .get("assembly", {})
                .get("molecule", [])
            )
            for molecule in (
                asu_molecules if isinstance(asu_molecules, list) else [asu_molecules]
            ):
                if not molecule:
                    continue

                # Ligand molecule
                if id_is_ligand(molecule["chain_id"]):
                    (auth_asym_id, ccd_id, auth_seq_id) = extract_ligand_contents(
                        molecule["chain_id"]
                    )

                    ligand_df = self._extract_ligand_rows(
                        ccd_id=ccd_id,
                        auth_asym_id=auth_asym_id,
                        auth_seq_id=auth_seq_id,
                    )

                    # Extract label_asym_id from dataframe
                    label_asym_id, label_seq_id = self._extract_ligand_ids(
                        ligand_df, auth_seq_id, auth_asym_id, ccd_id
                    )

                    molecule["label_seq_id_start"] = label_seq_id
                    molecule["label_seq_id_end"] = label_seq_id

                # Polymer molecule
                else:
                    auth_asym_id = molecule["chain_id"]
                    polymer_df = self._extract_polymer_rows(auth_asym_id)

                    # Extract label_asym_id from dataframe
                    (
                        label_asym_id,
                        auth_seq_id_start,
                        auth_seq_id_end,
                        label_seq_id_start,
                        label_seq_id_end,
                    ) = self._validate_polymer_ids(polymer_df, auth_asym_id)

                    molecule["label_asym_id"] = label_asym_id
                    molecule["auth_seq_id_start"] = auth_seq_id_start
                    molecule["auth_seq_id_end"] = auth_seq_id_end
                    molecule["label_seq_id_start"] = label_seq_id_start
                    molecule["label_seq_id_end"] = label_seq_id_end

        # Add missing fields for PDB files
        else:
            # ASM set (also known as PQS set)
            residue_mappings = get_residue_mappings(structure)

            for pqs_set in assembly_data["pisa_results"].get("asm_set", []):
                if not pqs_set:
                    continue

                if isinstance(pqs_set.get("assembly"), dict):
                    pqs_set["assembly"] = [pqs_set["assembly"]]

                for complex in pqs_set.get("assembly", []):
                    if not complex:
                        continue

                    for molecule in complex.get("molecule", []):
                        if not molecule:
                            continue

                        # Ligand molecule
                        if id_is_ligand(molecule["chain_id"]):
                            (
                                auth_asym_id,
                                ccd_id,
                                auth_seq_id,
                            ) = extract_ligand_contents(molecule["chain_id"])

                            if (
                                auth_asym_id in residue_mappings
                                and ccd_id in residue_mappings[auth_asym_id]["ligands"]
                            ):
                                molecule["auth_seq_id_start"] = auth_seq_id
                                molecule["auth_seq_id_end"] = auth_seq_id

                        # Polymer molecule
                        else:
                            auth_asym_id = molecule["chain_id"]
                            if auth_asym_id in residue_mappings:
                                molecule["auth_seq_id_start"] = residue_mappings[
                                    auth_asym_id
                                ]["polymers"]["start"]
                                molecule["auth_seq_id_end"] = residue_mappings[
                                    auth_asym_id
                                ]["polymers"]["end"]

            # ASU complex
            asu_molecules = (
                assembly_data["pisa_results"]
                .get("asu_complex", {})
                .get("assembly", {})
                .get("molecule", [])
            )

            for molecule in (
                asu_molecules if isinstance(asu_molecules, list) else [asu_molecules]
            ):
                if not molecule:
                    continue

                # Ligand molecule
                if id_is_ligand(molecule["chain_id"]):
                    (auth_asym_id, ccd_id, auth_seq_id) = extract_ligand_contents(
                        molecule["chain_id"]
                    )

                    if (
                        auth_asym_id in residue_mappings
                        and ccd_id in residue_mappings[auth_asym_id]["ligands"]
                    ):
                        molecule["auth_seq_id_start"] = auth_seq_id
                        molecule["auth_seq_id_end"] = auth_seq_id

                # Polymer molecule
                else:
                    auth_asym_id = molecule["chain_id"]
                    if auth_asym_id in residue_mappings:
                        molecule["auth_seq_id_start"] = residue_mappings[auth_asym_id][
                            "polymers"
                        ]["start"]
                        molecule["auth_seq_id_end"] = residue_mappings[auth_asym_id][
                            "polymers"
                        ]["end"]

        # Add CSS scores to interfaces in assembly
        for pqs_set in assembly_data["pisa_results"].get("asm_set", []):
            if not pqs_set:
                continue

            if isinstance(pqs_set.get("assembly"), dict):
                pqs_set["assembly"] = [pqs_set["assembly"]]

            for complex in pqs_set.get("assembly", []):
                if not complex:
                    continue

                interfaces = complex.get("interfaces", {}).get("interface", [])
                for interface in (
                    interfaces if isinstance(interfaces, list) else [interfaces]
                ):
                    interface_id = interface["id"]
                    interface_json_path = os.path.join(
                        self.path_interface_jsons,
                        f"interface_{interface_id}.json",
                    )

                    try:
                        with open(interface_json_path, "r") as interface_json_file:
                            interface_data = json.load(interface_json_file)
                            interface["css"] = interface_data.get("interface", {}).get(
                                "css", None
                            )
                    except FileNotFoundError:
                        LOGGER.warning(
                            f"Interface JSON file not found: {interface_json_path}. "
                            f"Skipping CSS addition."
                        )

        # Validate through pydantic model
        assembly_data = Complex(**assembly_data["pisa_results"]).model_dump()

        # Write JSON
        with open(self.path_output, "w") as json_file:
            # Parse via a pydantic model. Assemblies must be lists
            json.dump(assembly_data, json_file, indent=4)
            LOGGER.info(f"JSON file written successfully: {self.path_output}")


class ConvertListTextToJSON(ABC):
    def __init__(
        self, path_txt: str, path_json: str, pisa_data_type: str = None
    ) -> None:
        """
        Base class for converting PISA -list text files to JSON files.

        :param path_txt: Path to the input text file.
        :type path_txt: str
        :param path_json: Path to the output JSON file.
        :type path_json: str
        :param pisa_data_type: Type of PISA data (e.g., 'assemblies', 'interfaces',
            'monomers'), defaults to None
        :type pisa_data_type: str, optional
        """

        self.path_txt = path_txt

        if not path_json.endswith(".json"):
            path_json += f"/{pisa_data_type}_extended.json"

        self.path_json = path_json
        self.pisa_data_type = pisa_data_type.capitalize()

    def __str__(self):
        return (
            f"Convert{self.pisa_data_type}ListToJSON: "
            f"{self.path_txt} -> {self.path_json}"
        )

    def _clean_table_line(self, line: str) -> list[str]:
        """
        Extract relevant contents from a line of the table

        :param line: Line (row) of file to clean.
        :type line: str
        :return: List of cleaned contents from the line.
        :rtype: list[str]
        """
        return line.strip().replace("|", " ").split()

    def _end_of_table(self, line: str) -> bool:
        """
        Check if line is the end of the table. The formatting of the end of the table
        is unique from mid-table rules.

        :param line: Line (row) of file to check.
        :type line: str
        :return: True if line is end of table, False otherwise.
        :rtype: bool
        """
        return line.startswith("----") and "'" in line

    def _mid_table_rule(self, line: str) -> bool:
        """
        Check if line is a mid-table rule line. Mid-table rules are not the end of the
        table, but are used to separate merged rows.

        :param line: Line (row) of file to check.
        :type line: str
        :return: True if line is a mid-table rule, False otherwise.
        :rtype: bool
        """
        return line.startswith("----") and "+" in line

    def _find_start_of_table(self, lines: list[str]) -> int:
        """
        Find the index of the first line of the table in the text file.

        :param lines: Contents of the text file as a list of lines.
        :type lines: list[str]
        :raises Exception: If start of table cannot be found.
        :return: Index of the first line of the table.
        :rtype: int
        """

        start_table = None
        for i, line in enumerate(lines):
            if self._is_first_line_of_table(line):
                return i + 2

        if start_table is None:
            raise Exception(
                f"Could not find start of monomer table in file: {self.path_txt}"
            )

    def _is_first_line_of_table(self, line: str) -> bool:
        """
        Determine if the given line is the first line of the table. Only works if child
        process operates on single-table data. Can be overridden by subclasses if
        needed.

        :param line: Line of file to check.
        :type line: str
        :return: True if line is first line of table, False otherwise.
        :rtype: bool
        """
        return False

    @abstractmethod
    def parse(self) -> None:
        """
        Convert a text file with a list of items to a JSON file.
        """
        with open(self.path_txt, "r") as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]

        return lines


class ConvertInterfaceListToJSON(ConvertListTextToJSON):
    def __init__(self, path_txt: str, path_json: str) -> None:
        """
        Converts a PISA-generated interface -list text file into a single JSON file.

        :param path_txt: Path to the input interface data text file.
        :type path_txt: str
        :param path_json: Path to the output JSON file.
        :type path_json: str
        """
        super().__init__(path_txt, path_json, "interfaces")

    def _is_first_line_of_table(self, line: str) -> bool:
        return (
            "##" in line
            and "Id" in line
            and "Monomer1" in line
            and "Monomer2" in line
            and "Symmetry" in line
            and "Area" in line
            and "DeltaG" in line
            and "Nhb" in line
            and "Nsb" in line
            and "Nds" in line
        )

    def parse(self) -> None:
        """
        Convert interface -list text file to JSON file.
        """

        lines = super().parse()

        # Check for no interfaces
        if lines and lines[-1] == "NO INTERFACES FOUND":
            LOGGER.warning(
                f"No interfaces found in interface list file: {self.path_txt}. "
                "Not writing interfaces extended JSON file."
            )
            return

        # Parse table
        start_table = self._find_start_of_table(lines)
        interface_data = []
        for line in lines[start_table:]:
            if self._end_of_table(line):
                break
            parts = self._clean_table_line(line)

            if len(parts) < 10:
                LOGGER.warning(
                    f"Skipping invalid interface line: {line}. "
                    f"Expected at least 10 parts, got {len(parts)}"
                )
                continue

            interface = {
                "interface_id": parts[0],
                "serial_number": parts[1],
                "monomer_1": parts[2],
                "monomer_2": parts[3],
                "symmetry_operation": parts[4],
                "symmetry_id": parts[5],
                "area": parts[6],
                "delta_g": parts[7],
                "nhb": parts[8],
                "nsb": parts[9],
                "nds": parts[10],
            }
            interface_data.append(interface)

        interfaces_data = {"interfaces": interface_data}
        interfaces_data = InterfaceExtended(**interfaces_data).model_dump()

        with open(self.path_json, "w") as json_file:
            json.dump(interfaces_data, json_file, indent=4)
            LOGGER.info(f"Extended interface data JSON written: {self.path_json}")


class ConvertAssemblyListToJSON(ConvertListTextToJSON):
    def __init__(self, path_txt: str, path_json: str) -> None:
        """
        Converts a PISA-generated assembly -list text file into a single JSON file.

        :param path_txt: Path to the input assembly data text file.
        :type path_txt: str
        :param path_json: Path to the output JSON file.
        :type path_json: str
        """
        super().__init__(path_txt, path_json, "assemblies")

    def _is_first_line_of_first_table(self, line: str) -> bool:
        return (
            "Set" in line
            and "No" in line
            and "Size" in line
            and "Id" in line
            and "ASA" in line
            and "BSA" in line
            and "DGdiss0" in line
            and "mG0" in line
            and "Formula" in line
        )

    def _is_first_line_of_second_table(self, line: str) -> bool:
        return (
            "Size" in line
            and "Id" in line
            and "ASA" in line
            and "BSA" in line
            and "DGdiss0" in line
            and "mG0" in line
            and "Formula" in line
            and "Set" not in line
            and "No" not in line
        )

    def _row_continued(self, line: str) -> bool:
        """
        Check if the given line (row) is a continuation of the previous line (row) in
        the table.

        :param line: Line (row) of the file to check
        :type line: str
        :return: True if line is a continuation of the previous line, False otherwise
        :rtype: bool
        """

        contains_column_breaks = line.count("|") == 3
        contains_one_non_empty_column = len(self._clean_table_line(line)) == 1

        return contains_column_breaks and contains_one_non_empty_column

    def _extract_formula(self, parts: list[str], line_index: int) -> str | None:
        """
        Extract the formula from the given line and check if it continues on the next
        lines. If the line has 9 parts, the formula is in the 9th part.

        :param parts: Indidual columns of a row
        :type parts: list[str]
        :param line_index: Current line index in the text file, used to check next
            lines for continued formula
        :type line_index: int
        :return: The formula value from the table
        :rtype: str | None
        """
        formula = ""

        if len(parts) == 9:
            formula += parts[8]

            # Check next lines
            for next_line in self.lines[line_index + 1 :]:
                if self._mid_table_rule(next_line):
                    break

                if self._end_of_table(next_line):
                    break

                if self._row_continued(next_line):
                    next_parts = self._clean_table_line(next_line)
                    formula += next_parts[0]

        return formula if formula else None

    def parse(self) -> None:
        """
        Convert assembly -list text file to JSON file.
        """

        self.lines = super().parse()

        # Find indexes of tables
        starts_main_tables = []
        start_asm_table = None
        start_asu_table = None
        for i, line in enumerate(self.lines):
            if self._is_first_line_of_first_table(line):
                start_asm_table = i + 2
                starts_main_tables.append(start_asm_table)

            if self._is_first_line_of_second_table(line):
                start_asu_table = i + 2
                break

        if not start_asm_table:
            LOGGER.warning(
                "Could not find first ASM table in assembly list file: "
                f"{self.path_txt}. Not an issue if no assemblies were found or "
                "--as-is flag was used."
            )

        if not start_asu_table:
            raise Exception(
                f"Could not find ASU table in assembly list file: {self.path_txt}"
            )

        # Parse first table
        pqs_data = []
        for start_index in starts_main_tables:
            for i, line in enumerate(self.lines[start_index:]):
                if self._end_of_table(line):
                    break

                if self._mid_table_rule(line):
                    continue

                if self._row_continued(line):
                    continue

                parts = self._clean_table_line(line)

                if len(parts) < 8:
                    LOGGER.warning(
                        f"Skipping invalid PQS line: {line}. "
                        f"Expected at least 8 parts, got {len(parts)}"
                    )
                    continue

                formula = self._extract_formula(parts, start_index + i)

                pqs = {
                    "set": parts[0],
                    "number": parts[1],
                    "size": parts[2],
                    "id": parts[3],
                    "asa": parts[4],
                    "bsa": parts[5],
                    "dgdiss0": parts[6],
                    "mg0": parts[7],
                    "formula": formula,
                }
                pqs_data.append(pqs)

        # Parse second table
        asu_data = []
        for line in self.lines[start_asu_table:]:
            if self._end_of_table(line):
                break

            # Extract
            parts = self._clean_table_line(line)
            n_parts = len(parts)

            if n_parts == 7:
                formula = parts[6]

            elif n_parts < 6:
                LOGGER.warning(
                    f"Skipping invalid ASU line: {line}. "
                    f"Expected at least 6 columns, got {n_parts}"
                )
                continue

            elif n_parts == 6 and all(is_int_or_float(part) for part in parts):
                LOGGER.warning(f"Formula is missing for ASU line: {line}. ")
                formula = None

            else:
                LOGGER.error(
                    f"ASU line has unexpected format: {line}. "
                    f"Expected either 6 or 7 numeric columns, got {n_parts}"
                )
                raise Exception(
                    f"ASU line has unexpected format: {line}. "
                    f"Expected either 6 or 7 numeric columns, got {n_parts}"
                )

            asu = {
                "size": parts[0],
                "id": parts[1],
                "asa": parts[2],
                "bsa": parts[3],
                "dgdiss0": parts[4],
                "mg0": parts[5],
                "formula": formula,
            }
            asu_data.append(asu)

        extended_data = {
            "pqs_data": pqs_data,
            "asu_data": asu_data,
        }
        extended_data = ComplexExtended(**extended_data).model_dump()

        with open(self.path_json, "w") as json_file:
            json.dump(extended_data, json_file, indent=4)
            LOGGER.info(f"Extended assembly data JSON written: {self.path_json}")


class ConvertComponentsListToJSON(ConvertListTextToJSON):
    def __init__(self, path_txt: str, path_json: str) -> None:
        """
        Converts a PISA-generated monomer -list text file into a single JSON file.

        :param path_txt: Path to the input monomer data text file.
        :type path_txt: str
        :param path_json: Path to the output JSON file.
        :type path_json: str
        """
        super().__init__(path_txt, path_json, "monomers")

    def _is_first_line_of_table(self, line: str) -> bool:
        return (
            "##" in line
            and "Id" in line
            and "Monomer" in line
            and "Class" in line
            and "Nat" in line
            and "Nres" in line
            and "Sat" in line
            and "Sres" in line
            and "Area" in line
            and "DeltaG" in line
        )

    def parse(self) -> None:
        """
        Convert monomer -list text file to JSON file.
        """
        lines = super().parse()

        # Find start of table
        start_table = self._find_start_of_table(lines)

        # Parse table
        monomer_data = []
        for line in lines[start_table:]:
            if self._end_of_table(line):
                break
            parts = self._clean_table_line(line)

            if len(parts) < 9:
                LOGGER.warning(
                    f"Skipping invalid monomer line: {line}. "
                    f"Expected at least 8 parts, got {len(parts)}"
                )
                continue

            monomer = {
                "serial_number": parts[0],
                "monomer_id": parts[1],
                "chain_id": parts[2],
                "monomer_class": parts[3],
                "total_atoms": parts[4],
                "total_residues": parts[5],
                "surface_atoms": parts[6],
                "surface_residues": parts[7],
                "area": parts[8],
                "delta_g": parts[9] if len(parts) > 9 else None,
            }
            monomer_data.append(monomer)

        monomers_data = {"components": monomer_data}
        monomers_data = Components(**monomers_data).model_dump()

        with open(self.path_json, "w") as json_file:
            json.dump(monomers_data, json_file, indent=4)
            LOGGER.info(f"Extended monomer data JSON written: {self.path_json}")


class CompileInterfaceSummaryJSON:
    def __init__(
        self, path_interface_jsons: str, path_assembly_json: str, path_output_json: str
    ) -> None:
        """
        Reads individual interface JSON files and compiles a summary JSON.

        :param path_interface_jsons: Path to directory containing interface JSON files.
        :type path_interface_jsons: str
        :param path_assembly_json: Path to assembly JSON file.
        :type path_assembly_json: str
        :param path_output_json: Path to output summary JSON file.
        :type path_output_json: str
        """

        self.path_interface_jsons = path_interface_jsons
        self.path_assembly_json = path_assembly_json
        self.path_output_json = path_output_json

    def _load_assembly_json(self) -> dict:
        """
        Loads assembly JSON and creates a mapping of interface IDs to complex keys.

        :return: Mapping of interface IDs to complex keys.
        :rtype: dict
        """

        # Load assembly data
        assembly_to_interface_map = {}
        with open(self.path_assembly_json, "r") as f:
            d = json.load(f)
            pqs_sets = d.get("pqs_sets", [])

        if not pqs_sets:
            return {}

        for pqs_set in pqs_sets:
            for complex in pqs_set["complexes"]:
                complex_key = complex["complex_key"]
                interfaces = complex.get("interfaces", {})
                int_ids = [
                    interface["interface_id"]
                    for interface in interfaces.get("interfaces", [])
                ]
                for int_id in int_ids:
                    if int_id not in assembly_to_interface_map:
                        assembly_to_interface_map[int_id] = []
                    assembly_to_interface_map[int_id].append(complex_key)

        return assembly_to_interface_map

    def parse(self) -> None:
        """
        Converts individual interface JSON files into a summary JSON file.
        """

        interface_files = [
            f
            for f in os.listdir(self.path_interface_jsons)
            if f.startswith("interface_") and f.endswith(".json")
        ]
        if not interface_files:
            LOGGER.warning(
                "No interface JSON files found in directory: "
                f"{self.path_interface_jsons}. No summary JSON will be created."
            )
            return

        # Load assembly data
        assembly_to_interface_map = self._load_assembly_json()
        LOGGER.info(f"Loaded necessary assembly data from: {self.path_assembly_json}")

        # Prepare output data structure
        interface_summary = {}
        interface_summary["n_interfaces"] = len(interface_files)
        interface_summary["interface_types"] = []

        # Extract data
        int_type_tracker = set()
        for interface_file in interface_files:
            interface_json_path = os.path.join(
                self.path_interface_jsons, interface_file
            )

            with open(interface_json_path, "r") as f:
                d = json.load(f)

            interface_id = d["interface_id"]
            int_data = d["interface"]
            int_type = int_data["int_type"]
            mol1 = int_data["molecules"][0]
            mol2 = int_data["molecules"][1]

            summary_data = {
                "interface_id": interface_id,
                # NOTE: consider adding interface type here
                "auth_asym_id_1": mol1["auth_asym_id"],
                "int_natoms_1": mol1["int_natoms"],
                "int_nres_1": mol1["int_nres"],
                "auth_asym_id_2": mol2["auth_asym_id"],
                "int_natoms_2": mol2["int_natoms"],
                "int_nres_2": mol2["int_nres"],
                "int_area": int_data["int_area"],
                "int_solv_energy": int_data["int_solv_energy"],
                "pvalue": int_data["pvalue"],
                "css": int_data.get("css", None),
                # NOTE: Consider adding x_rel here
                "complex_keys_with_interface": assembly_to_interface_map.get(
                    interface_id, []
                ),
            }

            # Add to output dict
            if int_type not in int_type_tracker:
                interface_summary["interface_types"].append(
                    {"int_type": int_type, "interfaces": [summary_data]}
                )
                int_type_tracker.add(int_type)
            else:
                for int_type_entry in interface_summary["interface_types"]:
                    if int_type_entry["int_type"] == int_type:
                        int_type_entry["interfaces"].append(summary_data)
                        break

        # Order interface_types by int_type
        interface_summary["interface_types"] = sorted(
            interface_summary["interface_types"], key=lambda x: x["int_type"]
        )

        # Order the interfaces within each int_type by interface_id
        for int_type_entry in interface_summary["interface_types"]:
            int_type_entry["interfaces"] = sorted(
                int_type_entry["interfaces"],
                key=lambda x: int(x["interface_id"]),
            )

        interface_summary = InterfaceSummary(**interface_summary).model_dump()

        with open(self.path_output_json, "w") as json_file:
            json.dump(interface_summary, json_file, indent=4)
            LOGGER.info(f"Interface summary JSON written: {self.path_output_json}")
