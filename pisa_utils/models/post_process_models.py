from typing import Optional
from enum import Enum
from pydantic import BaseModel, Field, RootModel, field_validator

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
from pisa_utils.models.labels import (
    COMPLEXES_IN_PQS_SET_POST_PROC,
    INTERFACE_COMPONENT_P_VALUE,
    INTERFACE_COMPONENT_SOLVATION_ENERGY,
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
