from typing import Optional
from enum import Enum
from pydantic import Field, RootModel, field_validator, model_validator

from pisa_utils.constants import PRECISION_DP, STANDARD_DP
from pisa_utils.models.data_fields import (
    AuthAsymIdField,
    ComplexAccessibleSurfaceAreaField,
    ComplexBuriedSurfaceAreaField,
    ComplexCompositionField,
    ComplexCopiesInUnitCellField,
    ComplexDissociationEnergyField,
    ComplexEntropyChangeField,
    ComplexFormulaField,
    ComplexInterfaceEnergyField,
    ComplexKeyField,
    ComplexNumberMacromoleculesField,
    ComplexSignificanceScoreField,
    ComplexSymmetryNumberField,
    ComponentIsolatedSolvationEnergyField,
    ComponentNumSurfaceAtomsField,
    ComponentNumSurfaceResiduesField,
    ComponentTotalAtomsField,
    ComponentTotalResiduesField,
    ComponentTotalSurfaceAreaField,
    InterfaceAreaField,
    InterfaceIdField,
    InterfaceNumAtomsField,
    InterfaceNumResiduesField,
    InterfaceSolvationEnergyField,
    InterfaceTypeField,
    MoleculeClassField,
    PQSSetIdField,
    PValueField,
    SymmetryIdField,
    SymmetryOperationField,
    TotalInterfacesField,
)
from pisa_utils.models.data_models import StrictModel
from pisa_utils.models.labels import (
    AUTH_ASYM_ID,
    INTERFACE_AREA,
    INTERFACE_COMPONENT_P_VALUE,
    INTERFACE_COMPONENT_SOLVATION_ENERGY,
    INTERFACE_N_ATOMS,
    INTERFACE_N_RESIDUES,
    ATOM_LABEL,
    ATOM_LABEL_EXAMPLES,
    AUTH_SEQ_ID,
    BOND_DISTANCES,
    INTERFACE_BOND_TYPES,
    RESIDUE_3_LETTER_CODE,
    RESIDUE_3_LETTER_EXAMPLES,
    RESIDUE_ASA,
    RESIDUE_BSA,
    RESIDUE_SEQ_ID,
    RESIDUE_SOLVATION_ENERGY,
)


class PISAAnalysisType(str, Enum):
    PQS = "PQS set"
    ASU = "Asymmetric unit"


class ComplexTableRow(StrictModel):
    pqs_set_id: Optional[int] = PQSSetIdField(default=None)
    pisa_analysis_type: PISAAnalysisType = Field(
        default=PISAAnalysisType.PQS, description="Type of complex from PISA analysis"
    )
    complex_key: int = ComplexKeyField()
    formula: Optional[str] = ComplexFormulaField()
    composition: str = ComplexCompositionField()
    asa: float = ComplexAccessibleSurfaceAreaField()
    bsa: float = ComplexBuriedSurfaceAreaField()
    int_energy: float = ComplexInterfaceEnergyField()
    diss_energy: float = ComplexDissociationEnergyField()
    entropy: float = ComplexEntropyChangeField()
    mmsize: int = ComplexNumberMacromoleculesField()
    n_uc: int = ComplexCopiesInUnitCellField()
    symmetry_number: int = ComplexSymmetryNumberField()
    n_interfaces: Optional[int] = TotalInterfacesField(default=0)


class ComplexTable(RootModel[list[ComplexTableRow]]):
    pass


class InterfaceTableRow(StrictModel):
    int_type: int = InterfaceTypeField()
    interface_id: int = InterfaceIdField()
    auth_asym_id_1: str = AuthAsymIdField(
        description=f"{AUTH_ASYM_ID} for first molecule in interface",
    )
    int_natoms_1: int = InterfaceNumAtomsField(
        description=f"{INTERFACE_N_ATOMS} for first molecule in interface"
    )
    int_nres_1: int = InterfaceNumResiduesField(
        description=f"{INTERFACE_N_RESIDUES} for first molecule in interface"
    )
    auth_asym_id_2: str = AuthAsymIdField(
        description=f"{AUTH_ASYM_ID} for second molecule in interface",
    )
    int_natoms_2: int = InterfaceNumAtomsField(
        description=f"{INTERFACE_N_ATOMS} for second molecule in interface",
    )
    int_nres_2: int = InterfaceNumResiduesField(
        description=f"{INTERFACE_N_RESIDUES} for second molecule in interface",
    )
    int_area: float = Field(
        ..., description=INTERFACE_AREA, examples=[150.5, 300.75, 12.0]
    )
    int_solv_energy: float = InterfaceSolvationEnergyField()
    pvalue: float = PValueField()
    # Present when --as-is set to false
    css: Optional[float] = ComplexSignificanceScoreField()
    complex_keys_with_interface: Optional[list[int]] = Field(
        None,
        description="List of complex keys where this interface is present",
        examples=[1, 2, 3],
    )
    n_h_bonds: Optional[int] = Field(
        0,
        description="Number of hydrogen bonds in the interface",
        examples=[0, 2, 5],
    )
    n_salt_bridges: Optional[int] = Field(
        0,
        description="Number of salt bridges in the interface",
        examples=[0, 1, 3],
    )
    n_ss_bonds: Optional[int] = Field(
        0,
        description="Number of disulfide bonds in the interface",
        examples=[0, 1, 2],
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


class InterfaceTable(RootModel[list[InterfaceTableRow]]):
    pass


class InterfaceDetailsComponent(StrictModel):
    auth_asym_id: str = AuthAsymIdField()
    molecule_class: str = MoleculeClassField()
    symmetry_operation: str = SymmetryOperationField()
    symmetry_id: str = SymmetryIdField()

    int_natoms: int = InterfaceNumAtomsField()
    surface_natoms: int = ComponentNumSurfaceAtomsField()
    total_atoms: int = ComponentTotalAtomsField()

    int_nres: int = InterfaceNumResiduesField()
    surface_nres: int = ComponentNumSurfaceResiduesField()
    total_residues: int = ComponentTotalResiduesField()

    int_area: float = InterfaceAreaField()
    surface_area: float = ComponentTotalSurfaceAreaField()

    solv_energy: Optional[float] = ComponentIsolatedSolvationEnergyField(default=None)

    int_solv_energy: float = InterfaceSolvationEnergyField(
        description=INTERFACE_COMPONENT_SOLVATION_ENERGY
    )
    pvalue: float = PValueField(description=INTERFACE_COMPONENT_P_VALUE)

    @field_validator("int_area", "surface_area", "solv_energy")
    @classmethod
    def standard_round(cls, v):
        if v is None:
            return v
        return round(v, STANDARD_DP)

    @field_validator("int_solv_energy", "pvalue")
    @classmethod
    def prevision_round(cls, v):
        return round(v, PRECISION_DP)


class InterfaceDetails(StrictModel):
    interface_id: int = InterfaceIdField()
    int_type: int = InterfaceTypeField()
    css: Optional[float] = ComplexSignificanceScoreField()
    components: list[InterfaceDetailsComponent] = Field(
        description="List of data on the two components that form the interface",
    )

    @field_validator("css")
    @classmethod
    def prevision_round(cls, v):
        if v is None:
            return v
        return round(v, PRECISION_DP)


class ResidueRow(StrictModel):
    # For building table's primary key
    auth_asym_id: str = AuthAsymIdField()
    label_comp_id: str = Field(..., examples=["MET"], validation_alias="auth_comp_id")
    auth_seq_id: int = Field(..., description=AUTH_SEQ_ID, examples=[1, 2, 5])

    # Data
    bonds: Optional[str] = Field(
        None,
        description=INTERFACE_BOND_TYPES,
        examples=["H", "S", "D", "C", "HS", "HC", "HSDC", "CSHD"],
    )
    asa: float = Field(..., description=RESIDUE_ASA, examples=[158.71])
    bsa: float = Field(..., description=RESIDUE_BSA, examples=[123.57])
    pc_residue_buried: float = Field(
        0.0,
        description="Percentage of the residue that is buried by the interface only",
        examples=[77.86],
    )
    solv_energy: float = Field(
        ...,
        description=RESIDUE_SOLVATION_ENERGY,
        examples=[-1.23],
    )

    # For filtering
    is_interfacing_residue: bool = Field(
        ...,
        description="Indicates whether the residue is part of the interface",
    )
    is_solvent_accessible_residue: bool = Field(
        ...,
        description="Indicates whether the residue is solvent accessible",
    )
    is_inaccessible_residue: bool = Field(
        ...,
        description="Indicates whether the residue is inaccessible",
    )

    @model_validator(mode="after")
    def calculate_pc_residue_buried(self):
        """
        Calculate the percentage of the residue that is buried based on ASA and BSA.
        """
        if self.asa == 0:
            # Technically undefined, but we're setting to 0 for practical purposes
            self.pc_residue_buried = 0
        else:
            # Round to the nearest integer
            self.pc_residue_buried = round((self.bsa / self.asa) * 100, 0)
        return self

    @field_validator("asa", "bsa")
    @classmethod
    def standard_round(cls, v):
        return round(v, STANDARD_DP)

    @field_validator("solv_energy")
    @classmethod
    def precision_round(cls, v):
        return round(v, PRECISION_DP)

    @model_validator(mode="after")
    def determine_residue_accessibility(self):
        """
        Determine the accessibility status of the residue based on its properties.

        These bools are for table display purposes.
        """
        self.is_interfacing_residue = self.bsa > 0
        self.is_solvent_accessible_residue = self.asa > 0 and self.bsa == 0
        self.is_inaccessible_residue = self.asa == 0 and self.bsa == 0

        return self


class InterfaceComponentResidueTable(StrictModel):
    structure_1_residues: list[ResidueRow] = Field(
        ...,
        description=(
            "Structure 1's residue interface properties (not all residues will be at "
            "the interface)"
        ),
    )
    structure_2_residues: list[ResidueRow] = Field(
        ...,
        description=(
            "Structure 2's residues interface properties (not all residues will be at "
            "the interface)"
        ),
    )


class InterfacingBondRow(StrictModel):
    # First chain info
    auth_asym_id_1: str = AuthAsymIdField(
        description=f"{AUTH_ASYM_ID} for first molecule in bond",
    )
    label_comp_id_1: str = Field(
        ...,
        description=RESIDUE_3_LETTER_CODE,
        examples=RESIDUE_3_LETTER_EXAMPLES,
        validation_alias="auth_comp_id_1",
    )
    auth_seq_id_1: int = Field(..., description=RESIDUE_SEQ_ID, examples=["178"])
    auth_atom_id_1: str = Field(
        ...,
        description=ATOM_LABEL,
        examples=ATOM_LABEL_EXAMPLES,
    )

    # The actual metric
    dist: float = Field(
        ...,
        description=BOND_DISTANCES,
        examples=[2.8, 3.0, 3.5],
    )

    # Second chain info
    auth_asym_id_2: str = AuthAsymIdField(
        description=f"{AUTH_ASYM_ID} for second molecule in bond",
    )
    label_comp_id_2: str = Field(
        ...,
        description=RESIDUE_3_LETTER_CODE,
        examples=RESIDUE_3_LETTER_EXAMPLES,
        validation_alias="auth_comp_id_2",
    )
    auth_seq_id_2: int = Field(..., description=RESIDUE_SEQ_ID)
    auth_atom_id_2: str = Field(
        ...,
        description=ATOM_LABEL,
        examples=ATOM_LABEL_EXAMPLES,
    )

    @field_validator("dist")
    @classmethod
    def round_results(cls, v):
        return round(v, PRECISION_DP)


class InterfacingBondsTable(RootModel[list[InterfacingBondRow]]):
    """
    Object can be used for the hydrogen bond, salt bridge, disulfide bond, and covalent
    bond tables.

    Will need separate implementations if the bond type tables ever diverge in
    structure.
    """

    pass
