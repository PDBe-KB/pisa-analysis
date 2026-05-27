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
    COMPONENT_TOTAL_ATOMS,
    COMPONENT_TOTAL_RESIDUES,
    COMPOSITION,
    COPIES_IN_UNIT_CELL,
    FORMULA,
    INTERFACE_AREA,
    INTERFACE_CSS,
    INTERFACE_N_ATOMS,
    INTERFACE_N_RESIDUES,
    INTERFACE_NUMBER,
    INTERFACE_P_VALUE,
    INTERFACE_SOLVATION_ENERGY,
    INTERFACE_TOTAL,
    INTERFACE_TYPE,
    ISOLATED_COMPONENT_ASA,
    MOLECULE_CLASS,
    N_COMPONENT_SURFACE_ATOMS,
    N_COMPONENT_SURFACE_RESIDUES,
    NUM_MACROMOLECULES,
    PQS_SET_ID,
    SOLVATION_ENERGY_ISOLATED_STRUCTURE,
    SYMMETRY_ID,
    SYMMETRY_NUMBER,
    SYMMETRY_OPERATION,
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


def SymmetryOperationField(**kwargs) -> Field:
    defaults = {
        "description": SYMMETRY_OPERATION,
        "examples": ["x,y,z", "-X-1,Y,-Z+1/2", "X,Y,Z"],
    }
    return Field(**{**defaults, **kwargs})


def SymmetryIdField(**kwargs) -> Field:
    defaults = {
        "description": SYMMETRY_ID,
        "examples": ["0_555", "1_555", "3_455"],
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


def PValueField(**kwargs) -> Field:
    defaults = {
        "description": INTERFACE_P_VALUE,
        "examples": [0.01, 0.05, 0.1, 0.9],
    }
    return Field(**{**defaults, **kwargs})


def InterfaceNumAtomsField(**kwargs) -> Field:
    defaults = {
        "description": INTERFACE_N_ATOMS,
        "examples": [100, 200, 325],
    }
    return Field(**{**defaults, **kwargs})


def InterfaceNumResiduesField(**kwargs) -> Field:
    defaults = {
        "description": INTERFACE_N_RESIDUES,
        "examples": [10, 25, 50],
    }
    return Field(**{**defaults, **kwargs})


def InterfaceAreaField(**kwargs) -> Field:
    defaults = {
        "description": INTERFACE_AREA,
        "examples": [150.5, 300.75, 12.0],
    }
    return Field(**{**defaults, **kwargs})


def InterfaceSolvationEnergyField(**kwargs) -> Field:
    defaults = {
        "description": INTERFACE_SOLVATION_ENERGY,
        "examples": [-18.2, -5.5, -10.0, -2.3],
    }
    return Field(**{**defaults, **kwargs})


def ComponentTotalSurfaceAreaField(**kwargs) -> Field:
    defaults = {
        "description": ISOLATED_COMPONENT_ASA,
        "examples": [11112.6],
    }
    return Field(**{**defaults, **kwargs})


def ComponentIsolatedSolvationEnergyField(**kwargs) -> Field:
    defaults = {
        "description": SOLVATION_ENERGY_ISOLATED_STRUCTURE,
        "examples": [-220.0],
    }
    return Field(**{**defaults, **kwargs})


def ComponentTotalAtomsField(**kwargs) -> Field:
    defaults = {
        "description": COMPONENT_TOTAL_ATOMS,
        "examples": [1846],
    }
    return Field(**{**defaults, **kwargs})


def ComponentTotalResiduesField(**kwargs) -> Field:
    defaults = {
        "description": COMPONENT_TOTAL_RESIDUES,
        "examples": [248],
    }
    return Field(**{**defaults, **kwargs})


def ComponentNumSurfaceAtomsField(**kwargs) -> Field:
    defaults = {
        "description": N_COMPONENT_SURFACE_ATOMS,
        "examples": [500, 1000, 1500],
    }
    return Field(**{**defaults, **kwargs})


def ComponentNumSurfaceResiduesField(**kwargs) -> Field:
    defaults = {
        "description": N_COMPONENT_SURFACE_RESIDUES,
        "examples": [50, 100, 150],
    }
    return Field(**{**defaults, **kwargs})
