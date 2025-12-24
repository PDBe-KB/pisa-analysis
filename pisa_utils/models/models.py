from enum import Enum


class LigandProcessingMode(str, Enum):
    AUTO = "auto"
    FIXED = "fixed"
    FREE = "free"


class DataModes(str, Enum):
    ASSEMBLIES = "assemblies"
    INTERFACES = "interfaces"
    MONOMERS = "monomers"
