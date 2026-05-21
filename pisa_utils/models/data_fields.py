from pydantic import Field

from pisa_utils.models.labels import (
    AUTH_ASYM_ID,
    AUTH_ASYM_ID_EXAMPLES,
    COMPLEX_ASA,
    COMPLEX_BSA,
    COMPLEX_DISS_ENERGY,
    COMPLEX_ENTROPY_CHANGE,
    COMPLEX_INSTANCE_ID,
    COMPLEX_INTERFACE_ENERGY,
    COMPOSITION,
    COPIES_IN_UNIT_CELL,
    FORMULA,
    INTERFACE_CSS,
    INTERFACE_NUMBER,
    INTERFACE_TOTAL,
    INTERFACE_TYPE,
    MOLECULE_CLASS,
    NUM_MACROMOLECULES,
    PQS_SET_ID,
    SYMMETRY_NUMBER,
)


def AuthAsymIdField(**kwargs) -> Field:
    defaults = {
        "description": AUTH_ASYM_ID,
        "examples": AUTH_ASYM_ID_EXAMPLES,
    }
    return Field(**{**defaults, **kwargs})


def PQSSetIdField(**kwargs) -> Field:
    defaults = {
        "description": PQS_SET_ID,
        "examples": [1, 2, 3, 25],
    }
    return Field(**{**defaults, **kwargs})


def ComplexKeyField(**kwargs) -> Field:
    defaults = {
        "description": COMPLEX_INSTANCE_ID,
        "examples": [1, 2, 3],
    }
    return Field(**{**defaults, **kwargs})


def ComplexFormulaField(**kwargs) -> Field:
    defaults = {
        "default": None,
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


def ComplexNumberMacromoleculesField(**kwargs) -> Field:
    defaults = {
        "description": NUM_MACROMOLECULES,
        "examples": [2, 6, 24],
    }
    return Field(**{**defaults, **kwargs})


def ComplexCopiesInUnitCellField(**kwargs) -> Field:
    defaults = {
        "description": COPIES_IN_UNIT_CELL,
        "examples": [0, 1, 2, 3],
    }
    return Field(**{**defaults, **kwargs})


def ComplexSymmetryNumberField(**kwargs) -> Field:
    defaults = {
        "description": SYMMETRY_NUMBER,
        "examples": [1, 2, 4, 6],
    }
    return Field(**{**defaults, **kwargs})


def TotalInterfacesField(**kwargs) -> Field:
    defaults = {
        "description": INTERFACE_TOTAL,
        "examples": [0, 1, 5, 58],
    }
    return Field(**{**defaults, **kwargs})


def InterfaceIdField(**kwargs) -> Field:
    defaults = {
        "description": INTERFACE_NUMBER,
        "examples": [1, 2, 3, 10],
    }
    return Field(**{**defaults, **kwargs})


def InterfaceTypeField(**kwargs) -> Field:
    defaults = {
        "description": INTERFACE_TYPE,
        "examples": [1, 2, 5],
    }
    return Field(**{**defaults, **kwargs})


def ComplexSignificanceScoreField(**kwargs) -> Field:
    defaults = {
        "default": None,
        "description": INTERFACE_CSS,
        "examples": [1.0, 0.8, 0.5],
    }
    return Field(**{**defaults, **kwargs})


def MoleculeClassField(**kwargs) -> Field:
    defaults = {
        "description": MOLECULE_CLASS,
        "examples": ["Protein", "DNA", "RNA", "Ligand"],
    }
    return Field(**{**defaults, **kwargs})
