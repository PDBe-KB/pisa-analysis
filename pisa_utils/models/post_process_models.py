from pydantic import BaseModel, Field, RootModel

from pisa_utils.models.data_fields import (
    ComplexAccessibleSurfaceAreaField,
    ComplexBuriedSurfaceAreaField,
    ComplexCompositionField,
    ComplexDissociationEnergyField,
    ComplexEntropyChangeField,
    ComplexFormulaField,
    ComplexInterfaceEnergyField,
    PQSSetIdField,
)
from pisa_utils.models.labels import COMPLEXES_IN_PQS_SET_POST_PROC


class ComplexTableRow(BaseModel):
    complex_key: int = Field()
    formula: str = ComplexFormulaField()
    composition: str = ComplexCompositionField()
    asa: float = ComplexAccessibleSurfaceAreaField()
    bsa: float = ComplexBuriedSurfaceAreaField()
    int_energy: float = ComplexInterfaceEnergyField()
    diss_energy: float = ComplexDissociationEnergyField()
    entropy: float = ComplexEntropyChangeField()


class PQSSetRow(BaseModel):
    pqs_set_id: int = PQSSetIdField()
    complexes: list[ComplexTableRow] = Field(
        default=[], description=COMPLEXES_IN_PQS_SET_POST_PROC
    )


class ComplexTable(RootModel[list[PQSSetRow]]):
    pass
