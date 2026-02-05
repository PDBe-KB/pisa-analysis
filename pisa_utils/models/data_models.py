import logging
from enum import Enum
from typing import Literal, Optional, Union

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    RootModel,
    field_validator,
    model_validator,
)

from pisa_utils.constants import EXCLUDE_SIFTS_XREF_DB_FIELDS, PRECISION_DP, STANDARD_DP
from pisa_utils.models.labels import (
    ASU_COMPLEX,
    ATOM_LABEL,
    ATOM_LABEL_EXAMPLES,
    AUTH_ASYM_ID,
    AUTH_ASYM_ID_EXAMPLES,
    AUTH_SEQ_ID,
    BOND_DISTANCES,
    CELL_I,
    CELL_J,
    CELL_K,
    CHEMICAL_POTENTIAL_EN,
    COMPLEX_ASA,
    COMPLEX_BSA,
    COMPLEX_DISS_ENERGY,
    COMPLEX_ENTROPY_CHANGE,
    COMPLEX_INSTANCE_ID,
    COMPLEX_STABILITY,
    COMPLEX_STABLE,
    COMPLEX_STANARD_DISS_ENERGY,
    COMPLEX_STANDARD_ENTROPY_CHANGE,
    COMPLEX_TYPE,
    COMPLEXES_IN_ASU,
    COMPLEXES_IN_PQS_SET,
    COMPONENT_NUMBER,
    COMPONENT_TOTAL_ATOMS,
    COMPONENT_TOTAL_RESIDUES,
    COMPONENT_TYPE_ID,
    COMPOSITION,
    COPIES_IN_UNIT_CELL,
    DISS_AREA,
    FIXED_INTERFACE,
    FORMULA,
    INSERTION_CODE,
    INTERFACE_AREA,
    INTERFACE_BOND_TYPES,
    INTERFACE_CONTAINS_COVALENT_LINKAGE,
    INTERFACE_COVALENT_BONDS,
    INTERFACE_CRYSTALLOGRAPHIC_CONTACT,
    INTERFACE_CSS,
    INTERFACE_ENERGY,
    INTERFACE_H_BONDS,
    INTERFACE_N_ATOMS,
    INTERFACE_N_BONDS,
    INTERFACE_N_H_BONDS,
    INTERFACE_N_RESIDUES,
    INTERFACE_N_SALT_BRIDGES,
    INTERFACE_N_SS_BONDS,
    INTERFACE_NUMBER,
    INTERFACE_OTHER_BONDS,
    INTERFACE_P_VALUE,
    INTERFACE_SALT_BRIDGES,
    INTERFACE_SOLVATION_ENERGY,
    INTERFACE_SS_BONDS,
    INTERFACE_TOTAL,
    INTERFACE_TYPE,
    ISOLATED_COMPONENT_ASA,
    JOB_STATUS,
    LABEL_ASYM_ID,
    LABEL_SEQ_ID,
    MOLECULE_CLASS,
    MULTIMERIC_STATE,
    N_COMPONENT_SURFACE_ATOMS,
    N_COMPONENT_SURFACE_RESIDUES,
    N_DISS,
    N_PQS_SETS,
    NUM_COMPONENTS,
    NUM_MACROMOLECULES,
    PDB_ID,
    PISA_STABILISATION_ENERGY,
    PQS_SET_ID,
    PQS_SETS,
    RESIDUE_3_LETTER_CODE,
    RESIDUE_3_LETTER_EXAMPLES,
    RESIDUE_ASA,
    RESIDUE_BSA,
    RESIDUE_SEQ_ID,
    RESIDUE_SERIAL_NUMBER,
    RESIDUE_SOLVATION_ENERGY,
    RXX,
    RXY,
    RXZ,
    RYX,
    RYY,
    RYZ,
    RZX,
    RZY,
    RZZ,
    SESSION_NAME,
    SOLVATION_ENERGY_ISOLATED_STRUCTURE,
    STABILITY_DESCR,
    STATUS,
    STATUS_DESCRIPTION,
    STATUS_NOTE,
    SYMMETRY_ID,
    SYMMETRY_NUMBER,
    SYMMETRY_OPERATION,
    SYMMETRY_OPERATION_NUMBER,
    TX,
    TY,
    TZ,
    VISUAL_ID,
    XRAY_RELATED,
)
from pisa_utils.utils import extract_ligand_contents, id_is_ligand

LOGGER = logging.getLogger(__name__)


def convert_single_obj_to_list(v: Union[dict, list], dtype=dict) -> list:
    """Convert single dict to list if needed."""
    if isinstance(v, dtype):
        return [v]
    return v


def convert_yes_no_to_bool(v: Literal["Yes", "yes", "No", "no"]) -> bool:
    """
    Convert a "Yes"/"No" string to boolean.

    :param v: String value
    :type v: Literal["Yes", "yes", "No", "no"]
    :raises ValueError: _description_
    :return: Boolean value
    :rtype: bool
    """
    if v is None or isinstance(v, bool):
        return v
    if v.lower() == "yes":
        return True
    elif v.lower() == "no":
        return False
    else:
        raise ValueError(f"Invalid boolean value: {v}")


def remove_internal_whitespace(v: str) -> str:
    """
    Remove extra internal whitespace within a string.

    :param v: String with possible extra internal whitespace
    :type v: str
    :return: String with single internal spaces only
    :rtype: str
    """
    if v is None:
        return v
    return " ".join(v.split())


class LigandPosition(Enum):
    AUTO = "auto"
    FIXED = "fixed"
    FREE = "free"


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)


class Bond(StrictModel):
    # First chain info
    auth_asym_id_1: str = Field(
        ...,
        description=AUTH_ASYM_ID,
        examples=AUTH_ASYM_ID_EXAMPLES,
        validation_alias="chain-1",
    )
    label_asym_id_1: Optional[str] = Field(
        None, description=LABEL_ASYM_ID, examples=[], validation_alias="label_asym_id-1"
    )

    orig_label_asym_id_1: Optional[str] = Field(
        None, validation_alias="orig_label_asym_id-1"
    )
    pdbx_sifts_xref_db_num_1: Optional[str] = Field(
        None,
        validation_alias="pdbx_sifts_xref_db_num-1",
        exclude=EXCLUDE_SIFTS_XREF_DB_FIELDS,
    )
    pdbx_sifts_xref_db_name_1: Optional[str] = Field(
        None,
        validation_alias="pdbx_sifts_xref_db_name-1",
        exclude=EXCLUDE_SIFTS_XREF_DB_FIELDS,
    )
    pdbx_sifts_xref_db_acc_1: Optional[str] = Field(
        None,
        validation_alias="pdbx_sifts_xref_db_acc-1",
        exclude=EXCLUDE_SIFTS_XREF_DB_FIELDS,
    )
    auth_comp_id_1: str = Field(
        ...,
        description=RESIDUE_3_LETTER_CODE,
        examples=RESIDUE_3_LETTER_EXAMPLES,
        validation_alias="res-1",
    )
    auth_seq_id_1: int = Field(
        ..., description=RESIDUE_SEQ_ID, examples=["178"], validation_alias="seqnum-1"
    )
    label_seq_id_1: Optional[int] = Field(
        None,
        description=LABEL_SEQ_ID,
        examples=[1, 2, 180],
        validation_alias="label_seqnum-1",
    )
    # FIXME - find structures with insertion codes to check they're populated
    inscode_1: Optional[str] = Field(
        None, description=INSERTION_CODE, validation_alias="inscode-1"
    )
    auth_atom_id_1: str = Field(
        ...,
        description=ATOM_LABEL,
        examples=ATOM_LABEL_EXAMPLES,
        validation_alias="atname-1",
    )

    # Second chain info
    auth_asym_id_2: str = Field(
        ...,
        description=AUTH_ASYM_ID,
        examples=AUTH_ASYM_ID_EXAMPLES,
        validation_alias="chain-2",
    )
    label_asym_id_2: Optional[str] = Field(
        None, description=LABEL_ASYM_ID, validation_alias="label_asym_id-2"
    )
    orig_label_asym_id_2: Optional[str] = Field(
        None, validation_alias="orig_label_asym_id-2"
    )
    pdbx_sifts_xref_db_acc_2: Optional[str] = Field(
        None,
        validation_alias="pdbx_sifts_xref_db_acc-2",
        exclude=EXCLUDE_SIFTS_XREF_DB_FIELDS,
    )
    pdbx_sifts_xref_db_num_2: Optional[str] = Field(
        None,
        validation_alias="pdbx_sifts_xref_db_num-2",
        exclude=EXCLUDE_SIFTS_XREF_DB_FIELDS,
    )
    pdbx_sifts_xref_db_name_2: Optional[str] = Field(
        None,
        validation_alias="pdbx_sifts_xref_db_name-2",
        exclude=EXCLUDE_SIFTS_XREF_DB_FIELDS,
    )
    auth_comp_id_2: str = Field(
        ...,
        description=RESIDUE_3_LETTER_CODE,
        examples=RESIDUE_3_LETTER_EXAMPLES,
        validation_alias="res-2",
    )
    auth_seq_id_2: int = Field(
        ..., description=RESIDUE_SEQ_ID, validation_alias="seqnum-2"
    )
    label_seq_id_2: Optional[int] = Field(None, validation_alias="label_seqnum-2")
    inscode_2: Optional[str] = Field(
        None, description=INSERTION_CODE, validation_alias="inscode-2"
    )
    auth_atom_id_2: str = Field(
        ...,
        description=ATOM_LABEL,
        examples=ATOM_LABEL_EXAMPLES,
        validation_alias="atname-2",
    )
    dist: float = Field(
        ...,
        description=BOND_DISTANCES,
        examples=[2.8, 3.0, 3.5],
    )

    @field_validator("dist")
    @classmethod
    def round_results(cls, v):
        return round(v, PRECISION_DP)

    @field_validator("label_seq_id_1", "label_seq_id_2", mode="after")
    @classmethod
    def convert_overflowed_int_to_none(cls, v):
        if v == -2147483647:
            return None
        return v


class BondsInfo(StrictModel):
    n_bonds: int = Field(..., description=INTERFACE_N_BONDS, examples=[1, 2, 100])
    bonds: Optional[list[Bond]] = Field(
        [],
        description="List of extended bond information",
        validation_alias="bond",
        examples=[
            [
                Bond(
                    **{
                        "chain-1": "A",
                        "label_asym_id-1": "A",
                        "res-1": "ARG",
                        "seqnum-1": 45,
                        "label_seqnum-1": 45,
                        "atname-1": "NH1",
                        "chain-2": "B",
                        "label_asym_id-2": "B",
                        "res-2": "GLU",
                        "seqnum-2": 78,
                        "label_seqnum-2": 78,
                        "atname-2": "OE1",
                        "dist": 2.8,
                    }
                )
            ]
        ],
    )

    @field_validator("bonds", mode="before")
    @classmethod
    def ensure_bond_is_list(cls, v: Union[dict, list[dict]]) -> list[dict]:
        return convert_single_obj_to_list(v)


class Residue(StrictModel):
    residue_serial_number: int = Field(
        ...,
        description=RESIDUE_SERIAL_NUMBER,
        examples=[1, 2, 3, 100],
        validation_alias="ser_no",
    )
    auth_comp_id: str = Field(..., examples=["MET"], validation_alias="name")
    auth_seq_id: int = Field(
        ..., description=AUTH_SEQ_ID, examples=[1, 2, 5], validation_alias="seq_num"
    )
    label_seq_id: Optional[int] = Field(
        None, description=LABEL_SEQ_ID, validation_alias="label_seq_num"
    )
    ins_code: Optional[str] = Field(None, description=INSERTION_CODE)
    bonds: Optional[str] = Field(
        None,
        description=INTERFACE_BOND_TYPES,
        examples=["H", "S", "D", "C", "HS", "HC", "HSDC", "CSHD"],
    )
    asa: float = Field(..., description=RESIDUE_ASA, examples=[158.7123629])
    bsa: float = Field(..., description=RESIDUE_BSA, examples=[123.567])
    solv_energy: float = Field(
        ...,
        description=RESIDUE_SOLVATION_ENERGY,
        examples=[-1.23],
        validation_alias="solv_en",
    )

    @field_validator("asa", "bsa")
    @classmethod
    def standard_round(cls, v):
        return round(v, STANDARD_DP)

    @field_validator("solv_energy")
    @classmethod
    def precision_round(cls, v):
        return round(v, PRECISION_DP)

    @field_validator("label_seq_id", mode="after")
    @classmethod
    def convert_overflowed_int_to_none(cls, v):
        if v == -2147483647:
            return None
        return v


class Residues(StrictModel):
    residues: list[Residue] = Field(
        ...,
        description="List of residue information at interface",
        validation_alias="residue",
    )

    @field_validator("residues", mode="before")
    @classmethod
    def ensure_residues_is_list(cls, v: Union[dict, list[dict]]) -> list[dict]:
        return convert_single_obj_to_list(v)


class Molecule(StrictModel):
    mol_id: int = Field(..., validation_alias="id")
    auth_asym_id: str = Field(
        ...,
        description=AUTH_ASYM_ID,
        examples=AUTH_ASYM_ID_EXAMPLES,
        validation_alias="chain_id",
    )
    label_asym_id: Optional[str] = Field(None, description=LABEL_ASYM_ID, examples=[])
    ccd_id: Optional[str] = Field(None)

    auth_seq_id_start: int = Field(None)
    auth_seq_id_end: int = Field(None)

    label_seq_id_start: Optional[int] = Field(None)
    label_seq_id_end: Optional[int] = Field(None)

    molecule_class: str = Field(
        ..., description=MOLECULE_CLASS, examples=["Protein"], validation_alias="class"
    )
    symmetry_id: Optional[str] = Field(
        None, description=SYMMETRY_ID, examples=["0_555", "1_555"]
    )
    symmetry_operation_number: Optional[int] = Field(
        None,
        description=SYMMETRY_OPERATION_NUMBER,
        examples=[1, 2, 3],
        validation_alias="symop_no",
    )
    symmetry_operation: Optional[str] = Field(
        None,
        description=SYMMETRY_OPERATION,
        examples=["x,y,z"],
        validation_alias="symop",
    )
    cell_i: Optional[int] = Field(None, description=CELL_I, examples=[0, 1, 2])
    cell_j: Optional[int] = Field(None, description=CELL_J, examples=[0, 1, 2])
    cell_k: Optional[int] = Field(None, description=CELL_K, examples=[0, 1, 2])

    rxx: float = Field(..., description=RXX, examples=[1.0])
    rxy: float = Field(..., description=RXY, examples=[0.0])
    rxz: float = Field(..., description=RXZ, examples=[0.0])
    tx: float = Field(..., description=TX, examples=[0.0])
    ryx: float = Field(..., description=RYX, examples=[0.0])
    ryy: float = Field(..., description=RYY, examples=[1.0])
    ryz: float = Field(..., description=RYZ, examples=[0.0])
    ty: float = Field(..., description=TY, examples=[0.0])
    rzx: float = Field(..., description=RZX, examples=[0.0])
    rzy: float = Field(..., description=RZY, examples=[0.0])
    rzz: float = Field(..., description=RZZ, examples=[1.0])
    tz: float = Field(..., description=TZ, examples=[0.0])

    int_natoms: int = Field(
        ..., description=INTERFACE_N_ATOMS, examples=[100, 200, 325]
    )
    int_nres: int = Field(..., description=INTERFACE_N_RESIDUES, examples=[10, 25, 50])
    int_area: float = Field(
        ..., description=INTERFACE_AREA, examples=[150.5, 300.75, 12.0]
    )
    int_solv_energy: float = Field(
        ...,
        description=INTERFACE_SOLVATION_ENERGY,
        examples=[-5.5, -10.0, -2.3],
        validation_alias="int_solv_en",
    )
    pvalue: float = Field(
        ..., description=INTERFACE_P_VALUE, examples=[0.01, 0.05, 0.1, 0.9]
    )
    residues: Residues = Field(
        ...,
        description="List of residue information for residues at the interface",
        examples=[
            [
                Residue(
                    **{
                        "ser_no": 1,
                        "name": "MET",
                        "seq_num": 45,
                        "label_seq_num": 45,
                        "ins_code": None,
                        "bonds": "",
                        "asa": 158.7123629,
                        "bsa": 123.567,
                        "solv_en": -1.23,
                    }
                )
            ],
        ],
    )

    @field_validator(
        "int_area",
        "rxx",
        "rxy",
        "rxz",
        "tx",
        "ryx",
        "ryy",
        "ryz",
        "ty",
        "rzx",
        "rzy",
        "rzz",
        "tz",
    )
    @classmethod
    def standard_round(cls, v):
        return round(v, STANDARD_DP)

    @field_validator("int_solv_energy", "pvalue")
    @classmethod
    def prevision_round(cls, v):
        return round(v, PRECISION_DP)

    @model_validator(mode="after")
    def extract_ligand_fields(self):
        """
        Ensure that ligand fields are properly formatted.

        :return: Self with ligand fields populated if needed
        :rtype: Molecule
        """
        if self.ccd_id is not None:
            return self

        if not id_is_ligand(self.auth_asym_id):
            return self

        self.auth_asym_id, self.ccd_id, res_id = extract_ligand_contents(
            self.auth_asym_id
        )
        self.auth_seq_id_start, self.auth_seq_id_end = res_id, res_id

        return self

    @model_validator(mode="after")
    def add_symmetry_id_if_missing(self):
        """
        Add symmetry_id if missing based on symmetry_operation_number and cell indices.

        NOTE:
        This method assumes that symmetry_operation_number and cell_i, cell_j,
        cell_k are provided, symmetry_operation_number is the first part of symmetry_id
        from the MoleculeLabels model, and the cell_ indices are mapped to the second
        part of symmetry_id by adding 5 to each cell index (with None mapped to 5).
        This is not described in the PISA documentation but inferred from example data.

        :return: Self with symmetry_id populated if it was missing
        :rtype: Molecule
        """
        if self.symmetry_id is not None:
            return self

        if self.symmetry_operation_number is not None:
            self.symmetry_id = f"{self.symmetry_operation_number}_"
            for cell_n in (self.cell_i, self.cell_j, self.cell_k):
                if cell_n is not None:
                    self.symmetry_id += str(cell_n + 5)
                else:
                    self.symmetry_id += "5"

        return self


class InterfaceInfo(StrictModel):
    int_type: int = Field(..., description=INTERFACE_TYPE, validation_alias="type")
    n_occ: int = Field(...)
    int_area: float = Field(
        ..., description=INTERFACE_AREA, examples=[150.5, 300.75, 12.0]
    )
    int_solv_energy: float = Field(
        ...,
        description=INTERFACE_SOLVATION_ENERGY,
        examples=[-5.5, -10.0, -2.3],
        validation_alias="int_solv_en",
    )
    pvalue: float = Field(
        ..., description=INTERFACE_P_VALUE, examples=[0.01, 0.05, 0.1, 0.9]
    )
    stab_energy: float = Field(
        ..., description=PISA_STABILISATION_ENERGY, validation_alias="stab_en"
    )

    # Present when --as-is set to false
    css: Optional[float] = Field(
        None, description=INTERFACE_CSS, examples=[1.0, 0.8, 0.5]
    )
    overlap: Optional[str] = Field(None, description=None, examples=["No"])
    x_rel: Optional[bool] = Field(
        None, description=XRAY_RELATED, examples=[True, False], validation_alias="x-rel"
    )
    fixed: Optional[bool] = Field(
        None, description=FIXED_INTERFACE, examples=[True, False]
    )

    h_bonds: BondsInfo = Field(
        ..., description=INTERFACE_H_BONDS, examples=[], validation_alias="h-bonds"
    )
    salt_bridges: BondsInfo = Field(
        ...,
        description=INTERFACE_SALT_BRIDGES,
        examples=[],
        validation_alias="salt-bridges",
    )
    ss_bonds: BondsInfo = Field(
        ..., description=INTERFACE_SS_BONDS, examples=[], validation_alias="ss-bonds"
    )
    cov_bonds: BondsInfo = Field(
        ...,
        description=INTERFACE_COVALENT_BONDS,
        examples=[],
        validation_alias="cov-bonds",
    )
    other_bonds: Optional[BondsInfo] = Field(
        None,
        description=INTERFACE_OTHER_BONDS,
        examples=[],
        validation_alias="other-bonds",
    )

    molecules: list[Molecule] = Field(..., validation_alias="molecule")

    @field_validator("pvalue", "css")
    @classmethod
    def precision_round(cls, v):
        return round(v, PRECISION_DP)

    @field_validator("int_area", "int_solv_energy", "stab_energy")
    @classmethod
    def standard_round(cls, v):
        return round(v, STANDARD_DP)

    @field_validator("x_rel", "fixed", mode="before")
    @classmethod
    def convert_str_to_bool(cls, v):
        return convert_yes_no_to_bool(v)


class Interface(StrictModel):
    interface_id: int = Field(..., description=INTERFACE_NUMBER, examples=[1, 2, 3, 10])
    n_interfaces: int = Field(..., description=INTERFACE_TOTAL, examples=[58, 100, 200])

    status: str = Field(..., description=STATUS, examples=["Ok"])

    pdb_id: Optional[str] = Field(None, description=PDB_ID, examples=["1atp", "6gve"])
    interface: InterfaceInfo = Field(
        ...,
        description="Interface information from PISA job",
    )


class MoleculeLabels(StrictModel):
    auth_asym_id: str = Field(
        ...,
        description=AUTH_ASYM_ID,
        examples=AUTH_ASYM_ID_EXAMPLES,
        validation_alias="chain_id",
    )
    label_asym_id: Optional[str] = Field(None, description=LABEL_ASYM_ID, examples=[])
    visual_id: Optional[str] = Field(
        None, description=VISUAL_ID, examples=["A", "B", "C", "a", "b", "c", "{-}"]
    )
    auth_seq_id_start: Optional[int] = Field(None)
    auth_seq_id_end: Optional[int] = Field(None)
    label_seq_id_start: Optional[int] = Field(None)
    label_seq_id_end: Optional[int] = Field(None)

    ccd_id: Optional[str] = Field(None)

    rxx: float = Field(..., description=RXX, examples=[1.0])
    rxy: float = Field(..., description=RXY, examples=[0.0])
    rxz: float = Field(..., description=RXZ, examples=[0.0])
    tx: float = Field(..., description=TX, examples=[0.0])
    ryx: float = Field(..., description=RYX, examples=[0.0])
    ryy: float = Field(..., description=RYY, examples=[1.0])
    ryz: float = Field(..., description=RYZ, examples=[0.0])
    ty: float = Field(..., description=TY, examples=[0.0])
    rzx: float = Field(..., description=RZX, examples=[0.0])
    rzy: float = Field(..., description=RZY, examples=[0.0])
    rzz: float = Field(..., description=RZZ, examples=[1.0])
    tz: float = Field(..., description=TZ, examples=[0.0])

    rxx_f: float = Field(..., validation_alias="rxx-f", examples=[0.0])
    rxy_f: float = Field(..., validation_alias="rxy-f", examples=[0.0])
    rxz_f: float = Field(..., validation_alias="rxz-f", examples=[0.0])
    tx_f: float = Field(..., validation_alias="tx-f", examples=[0.0])
    ryx_f: float = Field(..., validation_alias="ryx-f", examples=[0.0])
    ryy_f: float = Field(..., validation_alias="ryy-f", examples=[0.0])
    ryz_f: float = Field(..., validation_alias="ryz-f", examples=[0.0])
    ty_f: float = Field(..., validation_alias="ty-f", examples=[0.0])
    rzx_f: float = Field(..., validation_alias="rzx-f", examples=[0.0])
    rzy_f: float = Field(..., validation_alias="rzy-f", examples=[0.0])
    rzz_f: float = Field(..., validation_alias="rzz-f", examples=[0.0])
    tz_f: float = Field(..., validation_alias="tz-f", examples=[0.0])
    symmetry_id: str = Field(
        ...,
        description=SYMMETRY_ID,
        examples=["0_555", "1_555"],
        validation_alias="symId",
    )

    @model_validator(mode="after")
    def extract_ligand_fields(self):
        """
        Ensure that ligand fields are properly formatted.

        :return: Self with ligand fields populated if needed
        :rtype: MoleculeLabels
        """
        if self.ccd_id is not None:
            return self

        if not id_is_ligand(self.auth_asym_id):
            return self

        self.auth_asym_id, self.ccd_id, res_id = extract_ligand_contents(
            self.auth_asym_id
        )
        self.auth_seq_id_start, self.auth_seq_id_end = res_id, res_id

        return self

    @field_validator(
        "rxx",
        "rxy",
        "rxz",
        "tx",
        "ryx",
        "ryy",
        "ryz",
        "ty",
        "rzx",
        "rzy",
        "rzz",
        "tz",
        "rxx_f",
        "rxy_f",
        "rxz_f",
        "tx_f",
        "ryx_f",
        "ryy_f",
        "ryz_f",
        "ty_f",
        "rzx_f",
        "rzy_f",
        "rzz_f",
        "tz_f",
    )
    @classmethod
    def round_results(cls, v):
        return round(v, STANDARD_DP)


class InterfaceLabel(StrictModel):
    interface_id: int = Field(
        ..., description=INTERFACE_NUMBER, examples=[1, 2, 60], validation_alias="id"
    )
    dissociates: bool = Field(..., description=None, examples=[True, False])
    css: Optional[float] = Field(
        None, description=INTERFACE_CSS, examples=[0.0, 1.0, 0.8, 0.5]
    )

    @field_validator("dissociates", mode="before")
    @classmethod
    def convert_str_to_bool(cls, v):
        return convert_yes_no_to_bool(v)


class InterfaceLabels(StrictModel):
    """
    A minimal list of interfaces.
    """

    n_interfaces: int = Field(..., description=INTERFACE_TOTAL, examples=[0, 1, 5, 58])
    interfaces: Optional[list[InterfaceLabel]] = Field(
        [],
        description="List of minimal interfaces",
        examples=[
            [
                InterfaceLabel(id=1, dissociates="No"),
                InterfaceLabel(id=2, dissociates="Yes"),
                InterfaceLabel(id=3, dissociates="No"),
            ]
        ],
        validation_alias="interface",
    )

    @field_validator("interfaces", mode="before")
    @classmethod
    def ensure_interfaces_is_list(cls, v: Union[dict, list[dict]]) -> list[dict]:
        return convert_single_obj_to_list(v)


class InterfaceSummaryInfo(StrictModel):
    interface_id: int = Field(..., description=INTERFACE_NUMBER, examples=[1, 2, 3])
    auth_asym_id_1: str = Field(
        ...,
        description=f"{AUTH_ASYM_ID} for first molecule in interface",
        examples=AUTH_ASYM_ID_EXAMPLES,
    )
    int_natoms_1: int = Field(
        ...,
        description=f"{INTERFACE_N_ATOMS} for first molecule in interface",
        examples=[100, 200, 325],
    )
    int_nres_1: int = Field(
        ...,
        description=f"{INTERFACE_N_RESIDUES} for first molecule in interface",
        examples=[10, 25, 50],
    )
    auth_asym_id_2: str = Field(
        ...,
        description=f"{AUTH_ASYM_ID} for second molecule in interface",
        examples=AUTH_ASYM_ID_EXAMPLES,
        validation_alias="chain_id",
    )
    int_natoms_2: int = Field(
        ...,
        description=f"{INTERFACE_N_ATOMS} for second molecule in interface",
        examples=[100, 200, 325],
    )
    int_nres_2: int = Field(
        ...,
        description=f"{INTERFACE_N_RESIDUES} for second molecule in interface",
        examples=[10, 25, 50],
    )
    int_area: float = Field(
        ..., description=INTERFACE_AREA, examples=[150.5, 300.75, 12.0]
    )
    int_solv_energy: float = Field(
        ...,
        description=INTERFACE_SOLVATION_ENERGY,
        examples=[-5.5, -10.0, -2.3],
        validation_alias="int_solv_en",
    )
    pvalue: float = Field(
        ..., description=INTERFACE_P_VALUE, examples=[0.01, 0.05, 0.1, 0.9]
    )
    # Present when --as-is set to false
    css: Optional[float] = Field(
        None, description=INTERFACE_CSS, examples=[1.0, 0.8, 0.5]
    )
    complex_keys_with_interface: Optional[list[int]] = Field(
        None,
        description="List of complex keys where this interface is present",
        examples=[1, 2, 3],
    )

    @field_validator("pvalue", "css")
    @classmethod
    def precision_round(cls, v):
        if v is None:
            return v
        return round(v, PRECISION_DP)

    @field_validator("int_area", "int_solv_energy")
    @classmethod
    def standard_round(cls, v):
        return round(v, STANDARD_DP)


class InterfaceTypeLabel(StrictModel):
    int_type: int = Field(..., description=INTERFACE_TYPE)
    interfaces: list[InterfaceSummaryInfo] = Field(
        ...,
        description="List of interface summaries for this interface type",
    )


class InterfaceSummary(StrictModel):
    n_interfaces: int = Field(..., description=INTERFACE_TOTAL, examples=[58, 100, 200])
    interface_types: list[InterfaceTypeLabel] = Field(
        ...,
        description="List of interface type summaries",
    )


class ComplexInfo(StrictModel):
    """Data for a given predicted complex"""

    complex_key: Optional[int] = Field(
        None,
        description=COMPLEX_INSTANCE_ID,
        examples=[1, 2, 3],
        validation_alias="serial_no",
    )
    complex_type: int = Field(
        ...,
        description=COMPLEX_TYPE,
        examples=[1, 2, 3],
        validation_alias="id",
    )
    size: int = Field(..., description=NUM_COMPONENTS, examples=[24, 32, 48])
    mmsize: int = Field(..., description=NUM_MACROMOLECULES, examples=[2, 6, 24])
    freesize: int = Field(..., description=None, examples=[16, 26, 40])
    stability_description: Optional[str] = Field(
        None,
        description=STABILITY_DESCR,
        examples=["This assembly appears to be stable in solution."],
        validation_alias="score",
    )
    diss_energy: float = Field(
        ...,
        description=COMPLEX_DISS_ENERGY,
        examples=[2.2527918698, 215.83686019],
    )
    diss_energy_0: Optional[float] = Field(
        None, description=COMPLEX_STANARD_DISS_ENERGY
    )
    stable: Optional[bool] = Field(
        None, description=COMPLEX_STABLE, examples=[True, False]
    )
    asa: float = Field(..., description=COMPLEX_ASA, examples=[154727.96462])
    bsa: float = Field(..., description=COMPLEX_BSA, examples=[64624.674488])
    entropy: float = Field(
        ..., description=COMPLEX_ENTROPY_CHANGE, examples=[76.165107939]
    )
    entropy_0: Optional[float] = Field(
        None, description=COMPLEX_STANDARD_ENTROPY_CHANGE
    )
    diss_area: float = Field(..., description=DISS_AREA, examples=[7375.0457086])
    int_energy: float = Field(
        ...,
        description=INTERFACE_ENERGY,
        examples=[-318.20296299],
    )
    n_uc: int = Field(
        ...,
        description=COPIES_IN_UNIT_CELL,
        examples=[0, 1, 2, 3],
    )
    n_diss: int = Field(
        ...,
        description=N_DISS,
        examples=[0, 1, 2, 3],
    )
    symmetry_number: int = Field(
        ..., description=SYMMETRY_NUMBER, examples=[4], validation_alias="symNumber"
    )
    formula: Optional[str] = Field(
        None, description=FORMULA, examples=[None, "A", "(2)", "A(8)B(4)C(4)a(8)"]
    )
    composition: str = Field(
        ...,
        description=COMPOSITION,
        examples=["ADEFIKNOBGJLCHMP[NAD](8)"],
    )
    interfaces: Optional[InterfaceLabels] = Field(
        None,
        description="Summary information of interfaces associated with the complex",
        examples=[
            InterfaceLabels(
                n_interfaces=58,
                interface=[
                    InterfaceLabel(id=1, dissociates="No"),
                    InterfaceLabel(id=2, dissociates="Yes"),
                    InterfaceLabel(id=3, dissociates="No"),
                ],
            )
        ],
    )
    molecules: list[MoleculeLabels] = Field(
        [],
        description="List of molecular components associated with the complex",
        validation_alias="molecule",
        examples=[
            [
                MoleculeLabels(
                    **{
                        "chain_id": "A",
                        "visual_id": None,
                        "rxx": 1,
                        "rxy": 0,
                        "rxz": 0,
                        "tx": 0,
                        "ryx": 0,
                        "ryy": 1,
                        "ryz": 0,
                        "ty": 0,
                        "rzx": 0,
                        "rzy": 0,
                        "rzz": 1,
                        "tz": 0,
                        "rxx-f": 1.0,
                        "rxy-f": 0.0,
                        "rxz-f": 0.0,
                        "tx-f": 0.0,
                        "ryx-f": 0.0,
                        "ryy-f": 1.0,
                        "ryz-f": 0.0,
                        "ty-f": 0.0,
                        "rzx-f": 0.0,
                        "rzy-f": 0.0,
                        "rzz-f": 1.0,
                        "tz-f": 0.0,
                        "symId": "1_555",
                    }
                ),
                MoleculeLabels(
                    **{
                        "chain_id": "B",
                        "visual_id": None,
                        "rxx": 1,
                        "rxy": 0,
                        "rxz": 0,
                        "tx": 0,
                        "ryx": 0,
                        "ryy": 1,
                        "ryz": 0,
                        "ty": 0,
                        "rzx": 0,
                        "rzy": 0,
                        "rzz": 1,
                        "tz": 0,
                        "rxx-f": 1.0,
                        "rxy-f": 0.0,
                        "rxz-f": 0.0,
                        "tx-f": 0.0,
                        "ryx-f": 0.0,
                        "ryy-f": 1.0,
                        "ryz-f": 0.0,
                        "ty-f": 0.0,
                        "rzx-f": 0.0,
                        "rzy-f": 0.0,
                        "rzz-f": 1.0,
                        "tz-f": 0.0,
                        "symId": "2_555",
                    }
                ),
            ],
        ],
    )

    @field_validator("molecules", mode="before")
    @classmethod
    def ensure_molecule_is_list(cls, v: Union[dict, list[dict]]) -> list[dict]:
        return convert_single_obj_to_list(v)

    @field_validator("diss_energy", "diss_energy_0")
    @classmethod
    def precision_round(cls, v):
        return round(v, PRECISION_DP)

    @field_validator("asa", "bsa", "entropy", "entropy_0", "diss_area", "int_energy")
    @classmethod
    def standard_round(cls, v):
        return round(v, STANDARD_DP)

    @field_validator("stability_description", mode="before")
    @classmethod
    def remove_extra_whitespace_within_string(cls, v):
        return remove_internal_whitespace(v)

    @model_validator(mode="after")
    def set_stable_flag(self):
        if self.diss_energy > 0:
            self.stable = True
        else:
            self.stable = False
        return self


class PQSSet(StrictModel):
    pqs_set_id: int = Field(
        ..., description=PQS_SET_ID, examples=[1, 2, 3, 25], validation_alias="ser_no"
    )
    all_chains_at_identity: bool = Field(..., description=None, examples=[True, False])
    stability: Optional[str] = Field(
        None,
    )
    complexes: list[ComplexInfo] = Field(
        ...,
        description=COMPLEXES_IN_PQS_SET,
        validation_alias="assembly",
        examples=[
            [
                ComplexInfo(
                    serial_no=1,
                    id=1,
                    size=2,
                    mmsize=2,
                    freesize=2,
                    score=("This assembly appears to be stable in solution."),
                    diss_energy=10.123456,
                    asa=12345.6789,
                    bsa=2345.6789,
                    entropy=50.123456,
                    diss_area=1234.5678,
                    int_energy=-150.12345,
                    n_uc=1,
                    n_diss=0,
                    symNumber=2,
                    formula="A(2)",
                    composition="AB(2)",
                    interfaces=None,
                    molecule=[],
                )
            ]
        ],
    )

    @field_validator("complexes", mode="before")
    @classmethod
    def ensure_complex_is_list(cls, v: Union[dict, list[dict]]) -> list[dict]:
        return convert_single_obj_to_list(v)

    @field_validator("all_chains_at_identity", mode="before")
    @classmethod
    def convert_yes_no_to_bool(cls, v):
        return convert_yes_no_to_bool(v)

    @model_validator(mode="after")
    def set_stability(self):
        # TODO the logic for determining stability is needed. Cannot find in docs.
        self.stability = None
        return self


class AsymmetricUnit(StrictModel):
    complex: ComplexInfo = Field(
        ...,
        description=COMPLEXES_IN_ASU,
        validation_alias="assembly",
        examples=[
            ComplexInfo(
                serial_no=1,
                id=1,
                size=2,
                mmsize=2,
                freesize=2,
                score=("This assembly appears to be stable in solution."),
                diss_energy=10.123456,
                asa=12345.6789,
                bsa=2345.6789,
                entropy=50.123456,
                diss_area=1234.5678,
                int_energy=-150.12345,
                n_uc=1,
                n_diss=0,
                symNumber=2,
                formula="A(2)",
                composition="AB(2)",
                interfaces=None,
                molecule=[],
            )
        ],
    )


class Complex(StrictModel):
    session_name: Optional[str] = Field(
        None,
        description=SESSION_NAME,
        examples=["1cbs", "e9ec45011a9da871da1a2c9cdc229723"],
        validation_alias="name",
    )

    status: str = Field(
        ...,
        description=JOB_STATUS,
        examples=["Ok"],
    )

    # Extant when --as-is is set
    status_description: Optional[str] = Field(
        None,
        description=STATUS_DESCRIPTION,
        examples=[
            (
                "Results are not available due to the following reason: Assemblies not "
                "calculated"
            )
        ],
    )
    status_note: Optional[str] = Field(
        None,
        description=STATUS_NOTE,
        examples=[
            (
                "If you think that this message is issued in error, please report it "
                "to author Eugene Krissinel at eugene.krissinel@stfc.ac.uk, providing "
                "all details of your input."
            )
        ],
    )

    # Extant when --as-is not set
    n_pqs_sets: Optional[int] = Field(
        None, description=N_PQS_SETS, examples=[1, 2, 3], validation_alias="total_asm"
    )
    n_interfaces: Optional[int] = Field(
        None, description=INTERFACE_TOTAL, examples=[0, 1, 5, 58]
    )
    assessment: Optional[str] = Field(
        None, description=COMPLEX_STABILITY, examples=["Stable", "Unstable"]
    )
    multimeric_state: Optional[int] = Field(
        None, description=MULTIMERIC_STATE, examples=[1, 2, 4, 8]
    )
    all_chains_at_identity: Optional[bool] = Field(
        None, description=None, examples=[True, False]
    )
    resolution: Optional[float] = Field(
        None,  # description=RESOLUTION, examples=[1.5, 2.0, 2.8]
    )
    pqs_sets: list[PQSSet] = Field(
        [],
        description=PQS_SETS,
        examples=[[PQSSet(ser_no=1, all_chains_at_identity="Yes", assembly=[])]],
        validation_alias="asm_set",
    )

    # Always extant
    asu_complex: Optional[AsymmetricUnit] = Field(
        None,
        description=ASU_COMPLEX,
        examples=[
            AsymmetricUnit(
                assembly=ComplexInfo(
                    serial_no=1,
                    id=1,
                    size=2,
                    mmsize=2,
                    freesize=2,
                    score=("This assembly appears to be stable in solution."),
                    diss_energy=10.123456,
                    asa=12345.6789,
                    bsa=2345.6789,
                    entropy=50.123456,
                    diss_area=1234.5678,
                    int_energy=-150.12345,
                    n_uc=1,
                    n_diss=0,
                    symNumber=2,
                    formula="A(2)",
                    composition="AB(2)",
                    interfaces=None,
                    molecule=[],
                )
            )
        ],
    )

    @field_validator("pqs_sets", mode="before")
    @classmethod
    def ensure_pqs_sets_is_list(cls, v: Union[dict, list[dict]]) -> list[dict]:
        return convert_single_obj_to_list(v)

    @field_validator("all_chains_at_identity", mode="before")
    @classmethod
    def convert_yes_no_to_bool(cls, v):
        return convert_yes_no_to_bool(v)

    @field_validator("status_note", mode="before")
    @classmethod
    def remove_extra_whitespace_within_string(cls, v):
        return remove_internal_whitespace(v)


class InterfaceExtensionLabels(StrictModel):
    interface_id: int = Field(..., description=INTERFACE_NUMBER, examples=[1, 2, 3])
    int_type: int = Field(
        ...,
        description=INTERFACE_TYPE,
        examples=[1, 2, 5],
        validation_alias="serial_number",
    )
    auth_asym_id_1: str = Field(
        ...,
        description=AUTH_ASYM_ID,
        examples=AUTH_ASYM_ID_EXAMPLES,
        validation_alias="monomer_1",
    )
    ccd_id_1: Optional[str] = Field(None)
    auth_seq_id_start_1: Optional[int] = Field(None)
    auth_seq_id_end_1: Optional[int] = Field(None)
    auth_asym_id_2: str = Field(
        ...,
        description=AUTH_ASYM_ID,
        examples=AUTH_ASYM_ID_EXAMPLES,
        validation_alias="monomer_2",
    )
    ccd_id_2: Optional[str] = Field(None)
    auth_seq_id_start_2: Optional[int] = Field(None)
    auth_seq_id_end_2: Optional[int] = Field(None)
    symmetry_operation: str = Field(
        ..., description=SYMMETRY_OPERATION, examples=["-X-1,Y,-Z+1/2"]
    )
    symmetry_id: str = Field(..., description=SYMMETRY_ID, examples=["3_455"])
    int_area: float = Field(
        ..., description=INTERFACE_AREA, examples=[1427.7], validation_alias="area"
    )
    int_solv_energy: float = Field(
        ...,
        description=INTERFACE_SOLVATION_ENERGY,
        examples=[-18.2],
        validation_alias="delta_g",
    )
    nhb: int = Field(..., description=INTERFACE_N_H_BONDS, examples=[0, 1, 2, 20])
    nsb: int = Field(..., description=INTERFACE_N_SALT_BRIDGES, examples=[0, 1, 2, 4])
    nds: int = Field(..., description=INTERFACE_N_SS_BONDS, examples=[0, 1, 2, 4])

    crystallographic_contact: Optional[bool] = Field(
        False, description=INTERFACE_CRYSTALLOGRAPHIC_CONTACT, examples=[True, False]
    )
    fixed_interface: Optional[bool] = Field(
        False, description=FIXED_INTERFACE, examples=[True, False]
    )
    contains_covalent_linkage: Optional[bool] = Field(
        False, description=INTERFACE_CONTAINS_COVALENT_LINKAGE, examples=[True, False]
    )

    @model_validator(mode="before")
    def extract_interface_property_from_id(self: dict):
        """
        Remove the trailing letter from interface ID if present. This letter was used
        to denote whether the interface was fixed (f) or crystal contact (x). However,
        this information is presented in the interface JSONs and is not needed here.

        NOTE: It's unclear what all the possible trailing characters are. From
        examples seen so far, these are:
         - f : fixed interface
         - x : crystal contact
         - c : interface contains covalent bonds
         - # : interface was fixed and also is a crystal contact (placeholder for "fx")
        """
        if isinstance(self["serial_number"], int):
            return self

        if self["serial_number"][-1].isdigit():
            return self

        characteristic = self["serial_number"][-1]
        match characteristic:
            case "f":
                self["fixed_interface"] = True
            case "x":
                self["crystallographic_contact"] = True
            case "c":
                self["contains_covalent_linkage"] = True
            case "#":
                self["fixed_interface"] = True
                self["crystallographic_contact"] = True
            case _:
                raise ValueError(
                    f"Unknown interface characteristic '{characteristic}' "
                    f"in interface ID '{self['serial_number']}'"
                )

        self["serial_number"] = self["serial_number"][:-1]

        return self

    @model_validator(mode="after")
    def extract_ligand_fields(self):
        """
        Ensure that ligand fields are properly formatted.

        :return: Self with ligand fields populated if needed
        :rtype: InterfaceExtensionLabels
        """

        # First ligand
        if self.ccd_id_1 is not None:
            return self
        if not id_is_ligand(self.auth_asym_id_1):
            return self
        self.auth_asym_id_1, self.ccd_id_1, res_id = extract_ligand_contents(
            self.auth_asym_id_1
        )
        self.auth_seq_id_start_1, self.auth_seq_id_end_1 = res_id, res_id

        # Second ligand
        if self.ccd_id_2 is not None:
            return self
        if not id_is_ligand(self.auth_asym_id_2):
            return self
        self.auth_asym_id_2, self.ccd_id_2, res_id = extract_ligand_contents(
            self.auth_asym_id_2
        )
        self.auth_seq_id_start_2, self.auth_seq_id_end_2 = res_id, res_id

        return self


class InterfaceExtended(StrictModel):
    interfaces: list[InterfaceExtensionLabels] = Field(
        ...,
        description="List of extended interface information",
        examples=[
            [
                InterfaceExtensionLabels(
                    serial_number=1,
                    interface_id="1",
                    monomer_1="A",
                    monomer_2="B",
                    symmetry_operation="-X-1,Y,-Z+1/2",
                    symmetry_id="3_455",
                    area=1427.7,
                    delta_g=-18.2,
                    nhb=2,
                    nsb=1,
                    nds=0,
                ),
                InterfaceExtensionLabels(
                    serial_number=2,
                    interface_id="2",
                    monomer_1="A",
                    monomer_2="C",
                    symmetry_operation="X,Y,Z",
                    symmetry_id="1_555",
                    area=850.3,
                    delta_g=-10.5,
                    nhb=0,
                    nsb=0,
                    nds=0,
                ),
            ]
        ],
    )

    @field_validator("interfaces", mode="before")
    @classmethod
    def ensure_interface_is_list(cls, v: Union[dict, list[dict]]) -> list[dict]:
        return convert_single_obj_to_list(v)


class PQSEntry(StrictModel):
    pqs_set_id: Optional[int] = Field(
        None, description=PQS_SET_ID, examples=[1, 2, 3], validation_alias="set"
    )
    complex_key: Optional[int] = Field(
        None,
        description=COMPLEX_INSTANCE_ID,
        examples=[1, 2, 3, 51],
        validation_alias="number",
    )
    complex_type: int = Field(
        ..., description=COMPLEX_TYPE, examples=[1, 2, 5], validation_alias="id"
    )
    mmsize: int = Field(
        ...,
        description=NUM_MACROMOLECULES,
        examples=[0, 1, 2, 25],
        validation_alias="size",
    )
    asa: float = Field(..., description=COMPLEX_ASA, examples=[19370.0])
    bsa: float = Field(..., description=COMPLEX_BSA, examples=[2855.3])
    diss_energy: float = Field(
        ...,
        description=COMPLEX_DISS_ENERGY,
        examples=[15.6],
        validation_alias="dgdiss0",
    )
    chem_energy: float = Field(
        ..., description=None, examples=[7.8], validation_alias="mg0"
    )
    formula: str = Field(..., description=FORMULA, examples=["A(2)", "A(2)a(2)b(2)"])

    @field_validator("asa", "bsa", "diss_energy", "chem_energy")
    @classmethod
    def round_results(cls, v):
        return round(v, STANDARD_DP)


class ASUEntry(StrictModel):
    complex_type: int = Field(
        ..., description=COMPLEX_TYPE, examples=[1, 2, 41], validation_alias="id"
    )
    mmsize: int = Field(
        ...,
        description=NUM_MACROMOLECULES,
        examples=[0, 1, 2, 25],
        validation_alias="size",
    )
    asa: float = Field(..., description=COMPLEX_ASA, examples=[11112.6])
    bsa: float = Field(..., description=COMPLEX_BSA, examples=[0.0])
    diss_energy: float = Field(
        ..., description=COMPLEX_DISS_ENERGY, examples=[0.0], validation_alias="dgdiss0"
    )
    chem_energy: float = Field(
        ..., description=CHEMICAL_POTENTIAL_EN, examples=[0.0], validation_alias="mg0"
    )
    formula: str = Field(..., description=FORMULA, examples=["A"])

    @field_validator("asa", "bsa", "diss_energy", "chem_energy")
    @classmethod
    def round_results(cls, v):
        return round(v, STANDARD_DP)


class ComplexExtended(StrictModel):
    pqs_data: list[PQSEntry] = Field(
        ...,
        description="List of additional PQS assembly information",
        examples=[
            [
                PQSEntry(
                    set=1,
                    number=1,
                    size=2,
                    id=1,
                    asa=19370.0,
                    bsa=2855.3,
                    dgdiss0=15.6,
                    mg0=7.8,
                    formula="A(2)",
                ),
            ]
        ],
    )
    asu_data: list[ASUEntry] = Field(
        ...,
        description="List of additional ASU complex information",
        examples=[
            [
                ASUEntry(
                    size=1,
                    id=1,
                    asa=11112.6,
                    bsa=0.0,
                    dgdiss0=0.0,
                    mg0=0.0,
                    formula="A",
                ),
            ]
        ],
    )

    @field_validator("pqs_data", "asu_data", mode="before")
    @classmethod
    def ensure_field_is_list(cls, v: Union[dict, list[dict]]) -> list[dict]:
        return convert_single_obj_to_list(v)


class Component(StrictModel):
    mol_id: int = Field(
        ...,
        description=COMPONENT_NUMBER,
        examples=[1],
        validation_alias="serial_number",
    )
    molecule_type_id: int = Field(
        ..., description=COMPONENT_TYPE_ID, examples=[1], validation_alias="monomer_id"
    )
    auth_asym_id: str = Field(
        ...,
        description=AUTH_ASYM_ID,
        examples=AUTH_ASYM_ID_EXAMPLES,
        validation_alias="chain_id",
    )
    ccd_id: Optional[str] = Field(None)
    auth_seq_id_start: Optional[int] = Field(None)
    auth_seq_id_end: Optional[int] = Field(None)
    molecule_class: str = Field(
        ...,
        description=MOLECULE_CLASS,
        examples=["Protein"],
        validation_alias="monomer_class",
    )
    total_atoms: int = Field(..., description=COMPONENT_TOTAL_ATOMS, examples=[1846])
    total_residues: int = Field(
        ..., description=COMPONENT_TOTAL_RESIDUES, examples=[248]
    )
    n_surface_atoms: int = Field(
        ...,
        description=N_COMPONENT_SURFACE_ATOMS,
        examples=[1017],
        validation_alias="surface_atoms",
    )
    n_surface_residues: int = Field(
        ...,
        description=N_COMPONENT_SURFACE_RESIDUES,
        examples=[223],
        validation_alias="surface_residues",
    )
    asa: float = Field(
        ...,
        description=ISOLATED_COMPONENT_ASA,
        examples=[11112.6],
        validation_alias="area",
    )
    solv_energy: Optional[float] = Field(
        None,
        description=SOLVATION_ENERGY_ISOLATED_STRUCTURE,
        examples=[-220.0],
        validation_alias="delta_g",
    )

    @field_validator("asa", "solv_energy")
    @classmethod
    def round_results(cls, v):
        if v is not None:
            return round(v, STANDARD_DP)
        return v

    @model_validator(mode="after")
    def extract_ligand_fields(self):
        """
        Ensure that ligand fields are properly formatted.

        :return: Self with ligand fields populated if needed
        :rtype: Component
        """
        if self.ccd_id is not None:
            return self

        if not id_is_ligand(self.auth_asym_id):
            return self

        self.auth_asym_id, self.ccd_id, res_id = extract_ligand_contents(
            self.auth_asym_id
        )
        self.auth_seq_id_start, self.auth_seq_id_end = res_id, res_id

        return self


class Components(StrictModel):
    components: list[Component] = Field(
        ...,
        description="List of component information",
        examples=[
            [
                Component(
                    serial_number=1,
                    monomer_id=1,
                    chain_id="A",
                    monomer_class="Protein",
                    total_atoms=1846,
                    total_residues=248,
                    surface_atoms=1017,
                    surface_residues=223,
                    area=11112.6,
                    delta_g=-220.0,
                ),
                Component(
                    serial_number=2,
                    monomer_id=2,
                    chain_id="B",
                    monomer_class="Protein",
                    total_atoms=1500,
                    total_residues=200,
                    surface_atoms=900,
                    surface_residues=180,
                    area=9500.3,
                    delta_g=-180.5,
                ),
            ]
        ],
    )

    @field_validator("components", mode="before")
    @classmethod
    def ensure_components_is_list(cls, v: Union[dict, list[dict]]) -> list[dict]:
        return convert_single_obj_to_list(v)


# For wrapping root element
class InterfaceOutput(RootModel[Interface]):
    root: Interface


class InterfaceSummaryOutput(RootModel[InterfaceSummary]):
    root: InterfaceSummary


class ComplexOutput(RootModel[Complex]):
    root: Complex


class InterfaceExtendedOutput(RootModel[InterfaceExtended]):
    root: InterfaceExtended


class ComplexExtendedOutput(RootModel[ComplexExtended]):
    root: ComplexExtended


class MonomersOutput(RootModel[Components]):
    root: Components
