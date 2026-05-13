from pydantic import Field

from pisa_utils.models.labels import (
    COMPLEX_ASA,
    COMPLEX_BSA,
    COMPLEX_DISS_ENERGY,
    COMPLEX_ENTROPY_CHANGE,
    COMPLEX_INTERFACE_ENERGY,
    COMPOSITION,
    FORMULA,
    PQS_SET_ID,
)


def PQSSetIdField(**kwargs) -> Field:
    defaults = {
        "description": PQS_SET_ID,
        "examples": [1, 2, 3, 25],
    }
    return Field(**{**defaults, **kwargs})


def ComplexFormulaField(**kwargs) -> Field:
    defaults = {
        "description": FORMULA,
        "examples": [None, "A", "(2)", "A(8)B(4)C(4)a(8)"],
    }
    return Field(**{**defaults, **kwargs})


def ComplexCompositionField(**kwargs) -> Field:
    defaults = {
        "description": COMPOSITION,
        "examples": ["ADEFIKNOBGJLCHMP[NAD](8)"],
    }
    return Field(**{**defaults, **kwargs})


def ComplexAccessibleSurfaceAreaField(**kwargs) -> Field:
    defaults = {
        "description": COMPLEX_ASA,
        "examples": [154727.96462],
    }
    return Field(**{**defaults, **kwargs})


def ComplexBuriedSurfaceAreaField(**kwargs) -> Field:
    defaults = {
        "description": COMPLEX_BSA,
        "examples": [0, 64624.67],
    }
    return Field(**{**defaults, **kwargs})


def ComplexInterfaceEnergyField(**kwargs) -> Field:
    defaults = {
        "description": COMPLEX_INTERFACE_ENERGY,
        "examples": [-0.5, 0, 5.6, -318.2],
    }
    return Field(**{**defaults, **kwargs})


def ComplexDissociationEnergyField(**kwargs) -> Field:
    defaults = {
        "description": COMPLEX_DISS_ENERGY,
        "examples": [2.25, 215.83],
    }
    return Field(**{**defaults, **kwargs})


def ComplexEntropyChangeField(**kwargs) -> Field:
    defaults = {
        "description": COMPLEX_ENTROPY_CHANGE,
        "examples": [76.17],
    }
    return Field(**{**defaults, **kwargs})
