# PISA Interface descriptions
PQS_SET_ID = """Probable quaternary structure (PQS) set is a group of different
complexes that possibly make up the crystal. PQS sets are ordered by PISA's estimated
probability of forming the complex inside the crystal."""

COMPLEX_STABILITY = """Labels whether a complex is likely (unstable) or unlikely
(stable) to dissociate in solution."""

COMPLEX_STABLE = """Indicates whether the complex is predicted to be stable in solution
based on its dissociation energy."""

COMPLEXES_IN_PQS_SET = """List of extended information about all complexes in the PQS
set. Each complex in the PQS set makes up a solution within the crystal."""

COMPLEXES_IN_ASU = """Extended information on the components of the crystallographic
unit cell taken as is, without prediction of complexes into PQS sets"""

JOB_STATUS = "Job success/failure status from the server"

SESSION_NAME = """Session name. Defaults to PDB ID for existing PDB entries or the job
ID for user-uploaded structures."""

STATUS_DESCRIPTION = "Description from PISA about the job status"

STATUS_NOTE = "Additional information from PISA about the job status"

N_PQS_SETS = """Total number of predicted assemblies. Returned only when submitted with
the 'complexes and interfaces' option selected."""

MULTIMERIC_STATE = """Number of macromolecular subunits in the most probable quaternary
structure (PQS) set."""

PQS_SETS = "List of extended information for all predicted complexes per PQS set"

ASU_COMPLEX = """Complex information for the asymmetric unit (ASU) of the crystal.
Always present, regardless of whether predicted complexes are found or if the "only
interfaces" option is selected."""

COMPLEX_INSTANCE_ID = """Unique integer number sequentially assigned to each complex
found across all PQS sets. Does not carry any biological meaning. Structurally
equivalent complexes across different PQS sets will have different complex instance
IDs."""

COMPLEX_TYPE = """An ID for structurally equivalent complexes (i.e. those with
equivalent components making the same set of interfaces)"""

NUM_COMPONENTS = "Number of macromolecular and ligand components in the complex"

NUM_MACROMOLECULES = "Number of macromolecular components only in the complex"

STABILITY_DESCR = """Short description of the stability of the complex in solution"""

COMPLEX_DISS_ENERGY = """The free energy of complex dissociation, in kcal/mol"""

COMPLEX_STANARD_DISS_ENERGY = """The standard free energy of complex dissociation at
RTP, in kcal/mol"""

COMPLEX_ASA = """The total solvent-accessible surface area, in A^2"""

COMPLEX_BSA = """The solvent-accessible surface area of components buried upon
formation of the given complex, in A^2"""

COMPLEX_ENTROPY_CHANGE = """The rigid-body entropy change at dissociation, at RTP in
kcal/mol"""

COMPLEX_STANDARD_ENTROPY_CHANGE = """The standard rigid-body entropy change at
dissociation, at RTP in kcal/mol"""

DISS_AREA = """Dissociation Interface Area (A^2)"""

INTERFACE_ENERGY = """The solvation free energy gain upon formation of the complex, in
kcal/mol"""

COPIES_IN_UNIT_CELL = """The number of given complexes found in a unit cell"""

N_DISS = """Number of dissociating parts"""

SYMMETRY_NUMBER = """Symmetry number indicates the number of different, but equivalent
orientations of the complex, which can be obtained by rotation"""

FORMULA = """Describes the components present within the predicted complex. Capital
letters denote macromolecular monomeric units (protein, DNA and RNA), lowercase letters
denote ligands."""

COMPOSITION = """Composition indicates monomeric units found in the predicted complex.
Chain IDs (auth_asym_id) and ligand CCDs from the uploaded file (or PDB entry) are
stated here. The “{-}” notation denotes chains missing an ID in the submitted file. """

INTERFACE_TYPE = """ID for interfaces made with equivalent sets of atoms from
equivalent components found across the entire crystal"""

INTERFACE_AREA = """The difference in total accessible surface areas of isolated and
interfacing structures, divided by two, in A^2"""

INTERFACE_SOLVATION_ENERGY = """The solvation free energy gain upon interface
formation, in kcal/mol"""

INTERFACE_P_VALUE = """The probability of getting a lower than observed ΔiG, when the
interface atoms are picked randomly from the protein surface"""

INTERFACE_CSS = """Complexation Significance Score (CSS) indicates the significance of
an interface towards complex formation"""

INTERFACE_H_BONDS = "Hydrogen bonds' extended information"
INTERFACE_SALT_BRIDGES = "Salt bridges' extended information"
INTERFACE_SS_BONDS = "Disulfide bonds' extended information"
INTERFACE_COVALENT_BONDS = "Covalent bonds' extended information"
INTERFACE_OTHER_BONDS = "Other bonds' extened information"

INTERFACE_N_BONDS = "Number of bonds of bond type at the interface"
INTERFACE_N_H_BONDS = "Number of hydrogen bonds at interface"
INTERFACE_N_SALT_BRIDGES = "Number of salt bridges at interface"
INTERFACE_N_SS_BONDS = "Number of disulfide bonds at interface"
INTERFACE_N_COVALENT_BONDS = "Number of covalent bonds at interface"
INTERFACE_N_OTHER_BONDS = """Number of other contacts within a distance of 4 A and not
classified as any of the other bonds"""

MOLECULE_CLASS = """Type of molecule"""

INTERFACE_N_ATOMS = "Number of atoms at interface"
INTERFACE_N_RESIDUES = "Number of residues at interface"

INTERFACE_TOTAL = "Total number of interfaces across all (predicted) complexes"
INTERFACE_NUMBER = """Unique interface identifier. Several unique interface identifiers
may correspond to the same interface type."""

RESIDUE_ASA = """Accessible surface area: The area in which the bulk solvent can
contact the residue, in A^2. Treated as dissociated if at the interface. """

RESIDUE_BSA = """Buried surface area: The surface area of the residue contributing to
the interface, in A^2"""

RESIDUE_SOLVATION_ENERGY = """Solvation energy of the corresponding residue, in
kcal/mol."""

RESIDUE_SERIAL_NUMBER = """Sequentially numbered identifier for the residue. Has no
correspondence to residue numbering in the input file."""

AUTH_SEQ_ID = """The residue sequence identifier as given in the uploaded file. For
.cif files, field is fetched from _atom_site.auth_seq_id. For .pdb files, fetched from
the residue sequence number column."""

LABEL_SEQ_ID = """The residue sequence identifier as given in the uploaded file. Comes
from _atom_site.label_seq_id in .cif files. For .pdb files, this is equivalent to the
auth_seq_id field."""

INSERTION_CODE = """The author-defined insertion code for the residue, if available."""

INTERFACE_BOND_TYPES = """Labels whether the interfacing bond is a hydrogen bond (H),
salt bridge (S), disulfide bond (D) or covalent bond (C)"""

SYMMETRY_ID = """The position of the component in the crystal, relative to the
coordinates given in the uploaded file (or PDB entry)."""

SYMMETRY_OPERATION = """The operation applied to a component to position it into the
corresponding complex"""

SYMMETRY_OPERATION_NUMBER = """The ID taken from the first number before the '_' in the
symmetry_id field."""

# NOTE These are presumed descriptions based interogating the data
CELL_I = """Corresponds to the first index of the three-index crystallographic
symmetry ID after the '_' in the symmetry_id field, minus 5."""
CELL_J = """Corresponds to the second index of the three-index crystallographic
symmetry ID after the '_' in the symmetry_id field, minus 5."""
CELL_K = """Corresponds to the third index of the three-index crystallographic
symmetry ID after the '_' in the symmetry_id field, minus 5."""

RXX = "Rotation matrix element at row 1, column 1 (X-axis component in X direction)"
RXY = "Rotation matrix element at row 1, column 2 (Y-axis component in X direction)"
RXZ = "Rotation matrix element at row 1, column 3 (Z-axis component in X direction)"
RYX = "Rotation matrix element at row 2, column 1 (X-axis component in Y direction)"
RYY = "Rotation matrix element at row 2, column 2 (Y-axis component in Y direction)"
RYZ = "Rotation matrix element at row 2, column 3 (Z-axis component in Y direction)"
RZX = "Rotation matrix element at row 3, column 1 (X-axis component in Z direction)"
RZY = "Rotation matrix element at row 3, column 2 (Y-axis component in Z direction)"
RZZ = "Rotation matrix element at row 3, column 3 (Z-axis component in Z direction)"

TX = "Translation component along the X-axis"
TY = "Translation component along the Y-axis"
TZ = "Translation component along the Z-axis"

XRAY_RELATED = """The interface is crystallographically related. When
true, the interface's monomers simultaneously form more than one given interface
in the crystal."""

FIXED_INTERFACE = """The interface was kept "fixed" during complex search."""

SOLVATION_ENERGY_ISOLATED_STRUCTURE = """Solvation energy of chain folding, in
kcal/mol"""

COMPONENT_TOTAL_ATOMS = """Total number of atoms in the component"""

COMPONENT_TOTAL_RESIDUES = """Total number of residues in the component"""

ISOLATED_COMPONENT_ASA = """The solvent-accessible surface area of a component,
in A^2"""

N_COMPONENT_SURFACE_ATOMS = """Number of solvent-accessible surface atoms in the
component"""

N_COMPONENT_SURFACE_RESIDUES = """Number of solvent-accessible surface residues in the
component"""

STATUS = "Status of analysis calculations"

PDB_ID = "PDB identifier (if available)"

AUTH_ASYM_ID = "Chain identifier (auth_asym_id)"
LABEL_ASYM_ID = """Chain identifier (label_asym_id) from mmCIF file. Not available in
PDB files."""
AUTH_ASYM_ID_EXAMPLES = ["A", "B", "aB", "12c"]

VISUAL_ID = """PISA-defined chain ID used in the 'formula' field. Chains with
sufficiently similar structural similarity are assigned the same visual ID."""

COMPONENT_NUMBER = """Unique ID corresponding to component instance in the (crystal)
structure"""

COMPONENT_TYPE_ID = """ID for components found within the (crystal) structure.
Components identified as structurally equivalent by PISA are assigned the ID."""

PISA_STABILISATION_ENERGY = "Stabilization energy (kcal/mol)"

CHEMICAL_POTENTIAL_EN = "Standard chemical potential energy, kcal/mol"

RESIDUE_3_LETTER_CODE = "Residue three-letter code"
RESIDUE_3_LETTER_EXAMPLES = ["THR", "ALA", "HIS"]
RESIDUE_SEQ_ID = "Residue sequence ID"
ATOM_LABEL = "Atom label"
ATOM_LABEL_EXAMPLES = ["CA", "N", "O", "CB", "OG1", "OD1"]


# PISA BondInfo descriptions
BOND_DISTANCES = "Bond distances (A)"
