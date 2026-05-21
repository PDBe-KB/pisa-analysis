from typing import Optional
from enum import Enum
from pydantic import BaseModel, Field, RootModel

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
    InterfaceIdField,
    InterfaceTypeField,
    MoleculeClassField,
    PQSSetIdField,
    TotalInterfacesField,
)
from pisa_utils.models.labels import COMPLEXES_IN_PQS_SET_POST_PROC


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
    pass
    auth_asym_id: str = AuthAsymIdField()
    molecule_class: str = MoleculeClassField()
    # symmetry_operation
    # symmetry_id
    # int_natoms
    # int_nres
    # int_area
    # int_natoms
    # int_nres
    # int_area
    # asa
    # solv_energy
    # int_solv_energy
    # pvalue


class InterfaceDetails(BaseModel):
    interface_id: int = InterfaceIdField()
    int_type: int = InterfaceTypeField()
    css: Optional[float] = ComplexSignificanceScoreField()
    components: list[InterfaceDetailsComponent] = Field()


class ComplexTable(RootModel[list[PQSSetRow]]):
    pass
