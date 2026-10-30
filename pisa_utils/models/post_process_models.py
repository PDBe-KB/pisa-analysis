from typing import Optional
from enum import Enum
from pydantic import BaseModel, Field, RootModel, field_validator, model_validator

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
    AUTH_SEQ_ID,
    COMPLEXES_IN_PQS_SET_POST_PROC,
    INSERTION_CODE,
    INTERFACE_BOND_TYPES,
    INTERFACE_COMPONENT_P_VALUE,
    INTERFACE_COMPONENT_SOLVATION_ENERGY,
    LABEL_SEQ_ID,
    RESIDUE_ASA,
    RESIDUE_BSA,
    RESIDUE_SOLVATION_ENERGY,
)


class PISAAnalysisType(str, Enum):
    PQS = "PQS set"
    ASU = "Asymmetric unit"


class ComplexTableRow(BaseModel):
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


class PQSSetRow(BaseModel):
    pqs_set_id: Optional[int] = PQSSetIdField(default=None)
    pisa_analysis_type: PISAAnalysisType = Field(
        default=PISAAnalysisType.PQS, description="Type of complex from PISA analysis"
    )
    complexes: list[ComplexTableRow] = Field(
        default=[], description=COMPLEXES_IN_PQS_SET_POST_PROC
    )


class InterfaceDetailsComponent(BaseModel):
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
    # FIXME: Change to component-specific name, avoid confusion with interface at large
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


class InterfaceDetails(BaseModel):
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


class ComplexTable(RootModel[list[PQSSetRow]]):
    pass


class InterfaceDetailsList(RootModel[list[InterfaceDetails]]):
    pass


class ResidueRow(StrictModel):
    # For building table's primary key
    auth_asym_id: str = AuthAsymIdField()
    label_comp_id: str = Field(..., examples=["MET"])
    auth_seq_id: int = Field(..., description=AUTH_SEQ_ID, examples=[1, 2, 5])
    label_seq_id: Optional[int] = Field(None, description=LABEL_SEQ_ID)
    ins_code: Optional[str] = Field(None, description=INSERTION_CODE)

    # Data
    bonds: Optional[str] = Field(
        None,
        description=INTERFACE_BOND_TYPES,
        examples=["H", "S", "D", "C", "HS", "HC", "HSDC", "CSHD"],
    )
    asa: float = Field(..., description=RESIDUE_ASA, examples=[158.71])
    bsa: float = Field(..., description=RESIDUE_BSA, examples=[123.57])
    pc_residue_buried: float = Field(
        ...,
        description="Percentage of the residue that is buried by the interface only",
        examples=[77.86],
    )
    solv_energy: float = Field(
        ...,
        description=RESIDUE_SOLVATION_ENERGY,
        examples=[-1.23],
        validation_alias="solv_en",
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
            self.pc_residue_buried = round((self.bsa / self.asa) * 100, None)
        return self

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
        description="Structure 1's residues at the interface",
    )
    structure_2_residues: list[ResidueRow] = Field(
        ...,
        description="Structure 2's residues at the interface",
    )
