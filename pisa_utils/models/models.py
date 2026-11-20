from enum import Enum


class LigandProcessingMode(str, Enum):
    AUTO = "auto"
    FIXED = "fixed"
    FREE = "free"


class DataModes(str, Enum):
    ASSEMBLIES = "assemblies"
    INTERFACES = "interfaces"
    MONOMERS = "monomers"


class OutputFormats(str, Enum):
    JSON = ".json"
    JSON_GZ = ".json.gz"
    XML = ".xml"
    XML_GZ = ".xml.gz"


class AllowedModelFileFormats(str, Enum):
    # Text-based formats
    CIF = ".cif"
    PDB = ".pdb"
    MMCIF = ".mmcif"
    ENT = ".ent"

    # Binary formats
    BCIF = ".bcif"

    # Compressed formats
    CIF_GZ = ".cif.gz"
    PDB_GZ = ".pdb.gz"
    MMCIF_GZ = ".mmcif.gz"
    ENT_GZ = ".ent.gz"
