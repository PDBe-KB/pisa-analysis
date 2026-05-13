from pydantic import BaseModel, Field

from pisa_utils.models.data_fields import (
    ComplexAccessibleSurfaceAreaField,
    ComplexBuriedSurfaceAreaField,
    ComplexCompositionField,
    ComplexDissociationEnergyField,
    ComplexFormulaField,
    ComplexInterfaceEnergyField,
    PQSSetIdField,
)
from pisa_utils.models.labels import COMPLEXES_IN_PQS_SET_POST_PROC


class ComplexTableRow(BaseModel):
    complex_key: int = Field()
    formula: str = ComplexFormulaField()
    composition: str = ComplexCompositionField()
    asa: str = ComplexAccessibleSurfaceAreaField()
    bsa: str = ComplexBuriedSurfaceAreaField()
    int_energy: float = ComplexInterfaceEnergyField()
    diss_energy: float = ComplexDissociationEnergyField()


class ComplexTable(BaseModel):
    pqs_set_id: int = PQSSetIdField()
    complexes: list[ComplexTableRow] = Field(
        default=[], description=COMPLEXES_IN_PQS_SET_POST_PROC
    )
