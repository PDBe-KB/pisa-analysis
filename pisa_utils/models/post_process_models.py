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
from pisa_utils.models.data_models import StrictModel
from pisa_utils.models.labels import (
    AUTH_ASYM_ID,
    COMPLEXES_IN_PQS_SET_POST_PROC,
    INTERFACE_AREA,
    INTERFACE_COMPONENT_P_VALUE,
    INTERFACE_COMPONENT_SOLVATION_ENERGY,
    INTERFACE_N_ATOMS,
    INTERFACE_N_RESIDUES,
    INTERFACE_TOTAL,
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


class PQSSetGroup(BaseModel):
    pqs_set_id: Optional[int] = PQSSetIdField(default=None)
    pisa_analysis_type: PISAAnalysisType = Field(
        default=PISAAnalysisType.PQS, description="Type of complex from PISA analysis"
    )
    complexes: list[ComplexTableRow] = Field(
        default=[], description=COMPLEXES_IN_PQS_SET_POST_PROC
    )


class ComplexTable(RootModel[list[PQSSetGroup]]):
    pass


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


class InterfaceDetailsList(RootModel[list[InterfaceDetails]]):
    pass


class InterfaceTableRow(StrictModel):
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
        validation_alias="chain_id",
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
    int_solv_energy: float = InterfaceSolvationEnergyField(
        validation_alias="int_solv_en",
    )
    pvalue: float = PValueField()
    # Present when --as-is set to false
    css: Optional[float] = ComplexSignificanceScoreField()
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


class InterfaceTableGroup(StrictModel):
    int_type: int = InterfaceTypeField()
    interfaces: list[InterfaceTableRow] = Field(
        ...,
        description="List of interface summaries for this interface type",
    )


class InterfaceTable(StrictModel):
    n_interfaces: int = Field(..., description=INTERFACE_TOTAL, examples=[58, 100, 200])
    interface_types: list[InterfaceTableGroup] = Field(
        ...,
        description="List of interface type summaries",
    )
