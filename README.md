# Assembly interfaces analysis (pisa-analysis)

[![Python >=3.10](https://img.shields.io/badge/python-%3E%3D3.10-blue.svg)](https://www.python.org/downloads/)

`pisa-analysis` is a Python library and command-line application for analysing
macromolecular interfaces and assemblies with PISA. It runs PISA, converts its XML
and text results into JSON, and prepares tables of complexes, interfaces, and
residues for downstream applications.

The installed distribution is named `pisa-analysis`; the Python import package is
`pisa_utils`. Two CLI commands provide different workflows:

| Command | Use it for |
| --- | --- |
| `pisa_analysis` | Analyse the interfaces in a supplied assembly and write two assembly-specific JSON files, optionally enriched with UniProt mappings. |
| `pisa_analysis_service` | Run the full PISA workflow, including assembly prediction, individual interface files, extended component data, and summary tables. Originally intended for PDBe services. |

Both workflows are supported. The service command can also analyse a supplied
structure without predicting assemblies by using `--as_is`.

## Installation

### Python package

Python **3.10 or later** is required. To install from this repository:

```shell
git clone https://github.com/PDBe-KB/pisa-analysis
cd pisa-analysis

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .

pisa_analysis --help
pisa_analysis_service --help
```

Installing the project installs its Python dependencies and registers both CLI
commands. Installing only `requirements.txt` installs the dependencies but does
not register the commands. For an editable installation, use
`python -m pip install -e .` instead.

The Python dependencies include Gemmi, pandas, Pydantic, lxml, xmltodict,
xmlschema, and jsonschema; their supported versions are defined in
[pyproject.toml](pyproject.toml). The exported runtime dependencies are recorded in
[requirements.txt](requirements.txt).

### PISA executable and setup data

Running either analysis CLI requires a separately installed PISA executable and
its setup data. Installing this Python package does not install PISA. See the
[PISA source and compilation instructions](https://github.com/PDBe-KB/pisa#compiling-and-running-pisa),
or use the [Docker image](#run-with-docker), which includes PISA.

Configure the executable and setup directory using absolute paths:

```shell
export PATH="/absolute/path/to/pisa/build:$PATH"
export PISA_SETUP_DIR="/absolute/path/to/pisa/setup"
```

You can also pass `--pisa_binary /absolute/path/to/pisa` and
`--pisa_setup_dir /absolute/path/to/setup` on each invocation. The executable
otherwise defaults to `pisa` on `PATH`; the setup directory falls back to
`PISA_SETUP_DIR`.

The setup directory must contain PISA's `srs/`, `molref/`, and `pisastore/` data and
a configuration template named
[`pisa_cfg_tmp`](https://github.com/PDBe-KB/pisa/blob/main/setup/pisa_cfg_tmp).
The wrappers create a temporary configuration from that template. The two
wrappers currently expect different session-directory placeholders:

| Workflow | Session-directory placeholder in `pisa_cfg_tmp` |
| --- | --- |
| `pisa_analysis` / `run_pisalite()` | `path_dataroot` |
| `pisa_analysis_service` / `run_pisa_service()` | `/data` |

For the service workflow, prepare a template with `path_dataroot` replaced by
`/data` if starting from the linked source template. Both wrappers replace
`path_to_setup` with the supplied setup directory. If you use both workflows,
keep their templates in separate setup directories with access to the same PISA
data. The implementation is in [file_io.py](pisa_utils/file_io.py) and
[run_pisa.py](pisa_utils/run_pisa.py).

Parsing existing PISA XML, text, or JSON through the Python library only requires
the Python dependencies; it does not require a PISA executable.

## Command-line usage

Run `pisa_analysis --help` or `pisa_analysis_service --help` to see the arguments
for the installed version. Both commands require these paths:

| Argument | Meaning |
| --- | --- |
| `-i`, `--input_cif` | Local input structure file. PISA accepts mmCIF or PDB; the option name is retained for both formats. |
| `-o`, `--output_json` | Directory for generated JSON files. |
| `--output_xml` | Directory for generated PISA XML and, for the service workflow, extended text results. |

Create output directories before running, particularly the JSON directory. Use a
separate output directory for each input structure to keep results together.
Both commands accept `-v` / `--verbose` for debug logging, as well as the PISA
executable and setup options described above.

### Analyse a supplied assembly

Use `pisa_analysis` when the input already contains the assembly you want to
analyse. `--pdb_id` and `--assembly_id` label the output; they do not download a
structure or select an assembly from an entry.

From a checkout, the bundled example can be run as follows after configuring
PISA:

```shell
mkdir -p output/6nxr/json output/6nxr/xml

pisa_analysis \
  --input_cif example_data/6nxr-assembly1.cif.gz \
  --pdb_id 6nxr \
  --assembly_id 1 \
  --input_updated_cif example_data/6nxr_updated.cif \
  --output_json output/6nxr/json \
  --output_xml output/6nxr/xml
```

The additional arguments for this command are:

| Argument | Meaning |
| --- | --- |
| `--pdb_id` | Required entry identifier used in filenames and JSON keys. |
| `--assembly_id` | Required assembly identifier used in filenames and JSON. |
| `--input_updated_cif` | Optional updated mmCIF file containing UniProt mappings. Adds UniProt accessions and residue numbering to interface contacts when available. |
| `--force` | Run PISA again even when both XML files already exist. |

This workflow runs PISA in `--as-is` mode, without predicting assemblies. It
writes `interfaces.xml` and `assembly.xml` into the XML directory. If both XML
files already exist, it reuses them unless `--force` is supplied; the reuse check
only tests their existence. JSON conversion still runs, and the CLI still
requires the PISA executable and setup directory to be configured.

The parser reads both XMLs, extracts interface contacts and assembly properties,
and optionally reads UniProt mappings from the updated mmCIF using Gemmi. It
writes these files into the JSON directory:

| File | Contents |
| --- | --- |
| `<pdb_id>-assembly<assembly_id>-interfaces.json` | Assembly properties and interface details, including molecules, residues, and bond contacts. |
| `<pdb_id>-assembly<assembly_id>.json` | Simplified assembly properties and the number of retained interfaces. |

For the example above, the names are `6nxr-assembly1-interfaces.json` and
`6nxr-assembly1.json`. Existing examples are available in
[example_data/json_output](example_data/json_output).

### Generate full PISA results

Use `pisa_analysis_service` for assembly prediction and the full set of parsed
outputs. A full crystallographic mmCIF is the intended input when predicting
assemblies. This command does not require `--pdb_id` or `--assembly_id`.

```shell
mkdir -p output/service/json/extended_data output/service/xml

pisa_analysis_service \
  --input_cif example_data/6nxr_updated.cif \
  --output_json output/service/json \
  --output_xml output/service/xml \
  --lig_proc_mode auto
```

Service-specific options are:

| Argument | Default | Meaning |
| --- | --- | --- |
| `--as_is` | Off | Calculate interfaces in the supplied structure without predicting assemblies. |
| `--lig_proc_mode {auto,fixed,free}` | `auto` | Choose PISA's ligand processing mode. |
| `--exclude_ligands CCD [CCD ...]` | None | Exclude ligands identified by CCD codes, for example `--exclude_ligands GOL SO4`. |
| `-z`, `--compress_output` | Off | Gzip the generated XML, extended text, and JSON files. |
| `--run_all_pisa_commands` | Enabled | Extract and parse assemblies, interfaces, and monomers. This is already enabled by default; the CLI currently has no option to disable it. |

The service always runs PISA. It accepts the shared `--force` flag, but that flag
does not change service behaviour. `--input_updated_cif` belongs to the assembly
command and is not a service option.

Create the JSON `extended_data/` directory as shown above; the current service
does not create that directory before writing its extended JSON files.

The service runs PISA's `-analyse`, `-xml`, and `-list` operations, then parses and
validates the results using the source structure and Pydantic models. The XML
output directory contains:

```text
assembly.xml
interfaces.xml
extended_data/
  assemblies_extended.txt
  interfaces_extended.txt
  monomers_extended.txt
```

The JSON output directory contains:

| Path | Contents |
| --- | --- |
| `interfaces/interface_<id>.json` | One interface per file, including bonds, both molecules, and their residue properties. |
| `assemblies.json` | Probable quaternary structure (PQS) sets and asymmetric-unit complex information. |
| `interface_summary.json` | Interfaces grouped by interface type, with links to complex keys. |
| `extended_data/assemblies_extended.json` | Parsed extended assembly data from PISA's `-list` output. |
| `extended_data/interfaces_extended.json` | Parsed extended interface data. |
| `extended_data/monomers_extended.json` | Parsed component/monomer data, including surface and total atom/residue counts. |
| `complex_table.json` | Flat list of complexes, with `pisa_analysis_type` and `pqs_set_id` on each row; ASU rows have a null `pqs_set_id`. |
| `interface_details.json` | Interface details joined to component/monomer data. |
| `interface_table.json` | Flat list of interfaces, with `int_type` and hydrogen-bond, salt-bridge, and disulfide-bond counts on each row, ordered by type and interface ID. |
| `residue_tables/interface_<id>_residue_table.json` | Per-interface residue tables for both components, including buried percentages and accessibility flags. |

With `--compress_output`, these paths have an additional `.gz` suffix. Interface
files and some derived outputs are omitted when there are no interfaces; the
full output listing assumes interfaces are present.

## Expected JSON formats

The two CLI workflows produce different JSON formats. For service output, the
validated structures and field descriptions are defined in
[data_models.py](pisa_utils/models/data_models.py) and
[post_process_models.py](pisa_utils/models/post_process_models.py). Concrete
examples are available in [tests/data/expected_output](tests/data/expected_output).

For assembly-specific output from `pisa_analysis`, the top-level key is the
supplied PDB ID. The bundled simplified assembly JSON looks like this:

```json
{
  "6nxr": {
    "assembly_id": "1",
    "pisa_version": "2.0",
    "assembly": {
      "id": "1",
      "size": "6",
      "interface_count": 1,
      "score": "",
      "macromolecular_size": "2",
      "dissociation_energy": 15.61,
      "accessible_surface_area": 19395.3,
      "buried_surface_area": 3514.17,
      "entropy": 12.98,
      "dissociation_area": 1427.5,
      "solvation_energy_gain": -35.28,
      "number_of_uc": "0",
      "number_of_dissociated_elements": "2",
      "symmetry_number": "2",
      "formula": "A(2)a(2)b(2)",
      "composition": "A-2A[NA](2)[GOL](2)",
      "R350": ""
    }
  }
}
```

Additional documentation for the assembly-interface format is available in the
[PISA-lite API documentation](https://pisalite.docs.apiary.io/#reference/0/pisaqualifierjson/interaction-interface-data-per-pdb-assembly-entry).
The newer service outputs use the models linked above.

## Python library usage

Import the functions and classes from their defining `pisa_utils` modules. The
parsers write files when `.parse()` is called. Use string paths: path handling
across the library does not uniformly support `pathlib.Path` objects. Create
parent output directories before parsing.

### Run PISA from Python

`run_pisa_service()` runs the external executable and returns paths to its XML
and optional extended text files. JSON conversion is a separate step.

```python
import os

from pisa_utils.run_pisa import run_pisa_service

assembly_xml, interfaces_xml, assemblies_txt, interfaces_txt, monomers_txt = (
    run_pisa_service(
        input_cif="example_data/6nxr_updated.cif",
        xml_output_dir="output/python/xml",
        pisa_binary="pisa",
        pisa_setup_dir=os.environ["PISA_SETUP_DIR"],
        extended_data=True,
    )
)
```

The function also accepts `asis`, `lig_proc_mode`, `exclude_ligands`, and
`compress_output`. With `extended_data=False`, it returns only the two XML paths.
For the assembly-only PISA execution wrapper, use `run_pisalite()` from the same
module. `AnalysePisa` in [analyze.py](pisa_utils/analyze.py) provides the
assembly-specific XML-to-JSON conversion used by `pisa_analysis`.

### Parse existing XML without running PISA

This example uses the bundled XML and matching structure fixtures, so it can be
run from a checkout without installing PISA:

```python
from pathlib import Path

from pisa_utils.parsers import (
    CompileInterfaceSummaryJSON,
    ConvertAssemblyXMLToJSON,
    ConvertInterfaceXMLToJSONs,
)

output = Path("output/parsed")
output.mkdir(parents=True, exist_ok=True)
interfaces = str(output / "interfaces")
structure = "tests/data/mock_data/af3_heterodimer_with_atp.cif"
xml_dir = "tests/data/mock_data/af3_heterodimer_with_atp"

ConvertInterfaceXMLToJSONs(
    path_xml=f"{xml_dir}/interfaces.xml",
    path_jsons=interfaces,
    path_structure_file=structure,
).parse()

ConvertAssemblyXMLToJSON(
    path_xml=f"{xml_dir}/assembly.xml",
    path_json=str(output / "assemblies.json"),
    path_structure_file=structure,
    path_interface_jsons=interfaces,
).parse()

CompileInterfaceSummaryJSON(
    path_interface_jsons=interfaces,
    path_assembly_json=str(output / "assemblies.json"),
    path_output_json=str(output / "interface_summary.json"),
).parse()
```

Parse interfaces before assemblies and compile the summary after both. For your
own data, supply the same structure file that was used to generate the XML so
chain and residue identifiers can be mapped correctly. PDB inputs are also
handled; mmCIF-specific label identifiers are unavailable for PDB files.

For gzip inputs, pass `compressed=True` and the `.gz` input paths. This setting
applies to both input and output. Most writers append `.gz` to output JSON names;
`CompileInterfaceSummaryJSON` expects its `path_output_json` to include `.gz`
when compression is enabled.

### Build tables from parsed JSON

Continue from the preceding parsing example:

```python
from pisa_utils.post_process_parsers import (
    PostProcessComplexTable,
    PostProcessInterfaceComponentResidueTable,
    PostProcessInterfaceTable,
)

PostProcessComplexTable(
    input_json_path=str(output / "assemblies.json"),
    output_json_path=str(output / "complex_table.json"),
).parse()

PostProcessInterfaceTable(
    input_summary_interfaces_path=str(output / "interface_summary.json"),
    path_interface_jsons=interfaces,
    output_json_path=str(output / "interface_table.json"),
).parse()

PostProcessInterfaceComponentResidueTable(
    path_interfaces=interfaces,
    output_json_dir=str(output / "residue_tables"),
).parse()
```

The residue parser reads `interface_<id>.json` files from an `interfaces/`
directory and writes one `interface_<id>_residue_table.json` per interface. Each
file follows `InterfaceComponentResidueTable`, with `structure_1_residues` and
`structure_2_residues`, including buried percentages and accessibility flags.
It preserves component and residue order and includes residues outside the
interface. Pass `compressed=True` for gzip input and output. The
`pisa_analysis_service` CLI runs this parser automatically and creates the
`residue_tables/` output directory when interface files are present.

Additional parsers can be used independently:

| Parser | Inputs and purpose |
| --- | --- |
| `ConvertAssemblyListToJSON` | Converts PISA's extended assembly text using `path_txt` and `path_json`. |
| `ConvertInterfaceListToJSON` | Converts extended interface text using `path_txt` and `path_json`. |
| `ConvertComponentsListToJSON` | Converts extended monomer text using `path_txt` and `path_json`. |
| `PostProcessInterfaceDetailsList` | Joins `path_interfaces` to `path_monomers_json` and writes combined details to `output_json_path`. |

The list-text parsers are in `pisa_utils.parsers`; the details parser is in
`pisa_utils.post_process_parsers`. For the list parsers, `path_json` may be an
existing output directory or a filename ending in `.json`. With compression,
pass the `.json` filename and `compressed=True`; the writer appends `.gz`.

### Overlay covariation data

The optional covariation module reads the assembly-interface JSON produced by
`pisa_analysis` and matches contacts using UniProt accessions and residue numbers:

```shell
mkdir -p output/covariation

python -m pisa_utils.covariations_int \
  --input_json tests/data/input_json.json \
  --input_cov tests/data/F5HCH8_cov.csv \
  --acc_unp_ids F5HCH8 \
  --output_csv output/covariation
```

When matching pairs are found, it writes
`<pdb_id>_assembly<assembly_id>-interfaces_cov.csv`. The covariation CSV uses
`uniprot_residue_index_a`, `uniprot_residue_index_b`, `covariation_probability`,
and `covariation_score` columns. Matching is implemented for a single UniProt
accession; heteromeric matching across multiple accessions is not implemented.

## Run with Docker

Build the Docker image from this checkout; it includes PISA and its setup data:

```shell
docker build . -t pisa-analysis
```

The [Dockerfile](Dockerfile) installs the package, making both `pisa_analysis`
and `pisa_analysis_service` available as CLI commands inside the container.

From a checkout, run `pisa_analysis_service` with a writable bind mount:

```shell
mkdir -p output/docker/json/extended_data output/docker/xml

docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" \
  pisa-analysis \
  pisa_analysis_service \
  --input_cif /data/example_data/6nxr_updated.cif \
  --output_json /data/output/docker/json \
  --output_xml /data/output/docker/xml \
  --lig_proc_mode auto
```

To run the assembly-specific workflow instead:

```shell
mkdir -p output/docker-assembly/json output/docker-assembly/xml

docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" \
  pisa-analysis \
  pisa_analysis \
  --input_cif /data/example_data/6nxr-assembly1.cif.gz \
  --pdb_id 6nxr \
  --assembly_id 1 \
  --input_updated_cif /data/example_data/6nxr_updated.cif \
  --output_json /data/output/docker-assembly/json \
  --output_xml /data/output/docker-assembly/xml
```

Replace `$PWD` with an absolute host directory if your files are elsewhere, and
use the corresponding `/data/...` paths inside the container. The examples use
your host UID and GID so generated files are owned by your user. Without
`--user`, the image runs as the non-root `app` user (UID 10001), which needs write
access to the mounted directories.

The [PISA base-image template](https://github.com/PDBe-KB/pisa/blob/v2.1.2/Dockerfile)
uses `/data` for session storage. The service substitutes a temporary session
directory; the assembly workflow keeps `/data`, so the mounted directory itself
must be writable and may also contain PISA session files.

## Development

The library supports Python 3.10 and later. The full development environment
requires Python 3.11 or later because pre-commit's transitive `filelock`
dependency requires it. CI runs the tests on Python 3.10, 3.11, and 3.12.

We use Astral's [`uv` tool](https://docs.astral.sh/uv/) to manage the project and
dependencies. From a checkout:

```shell
curl -LsSf https://astral.sh/uv/install.sh | sh
uv sync --locked --python 3.11
source .venv/bin/activate

pre-commit install
pre-commit run --all-files
make test
```

`uv sync` installs the project, its CLI commands, and the default development
dependency group, including the test tools. The pre-commit hooks keep
`requirements.txt` and `requirements-dev.txt` up to date and lint/format the code
with [Ruff](https://docs.astral.sh/ruff/). Dependency groups and tool configuration
live in [pyproject.toml](pyproject.toml) and
[.pre-commit-config.yaml](.pre-commit-config.yaml).

To run tests on Python 3.10 without the full development tools:

```shell
uv run --python 3.10 --no-dev --group test make test
```

For a direct or focused test run:

```shell
uv run --python 3.11 pytest
uv run --python 3.11 pytest tests/parsers/test_post_process_parsers.py
```

`make test` runs pytest with coverage and writes `report.xml` and
`coverage/cobertura-coverage.xml`. The tests use bundled XML, text, JSON, and
structure fixtures, so they do not need to execute PISA.

### Repository layout

| Path | Purpose |
| --- | --- |
| `pisa_utils/run.py`, `cli_tools.py` | CLI entry points and argument handling. |
| `pisa_utils/run_pisa.py` | PISA execution and extraction of XML/extended text. |
| `pisa_utils/analyze.py`, `dictionaries.py`, `field_handlers.py` | Assembly-specific JSON conversion, contact extraction, and identifier mappings. |
| `pisa_utils/parsers.py` | XML/text-to-JSON parsers and interface summaries. |
| `pisa_utils/post_process_parsers.py` | Derived complex, interface, and residue tables. |
| `pisa_utils/models/` | Pydantic models, field descriptions, and validation errors. |
| `pisa_utils/file_io.py` | JSON/gzip helpers and PISA configuration generation. |
| `pisa_utils/covariations_int.py` | Covariation matching and CSV output. |
| `example_data/` | Example structure, XML, and assembly-specific JSON files. |
| `tests/` | Unit tests and input/expected-output fixtures. |

## Troubleshooting

- **Command not found:** activate the virtual environment and install the project
  with `python -m pip install .`; installing dependencies alone does not install
  CLI commands.
- **PISA executable or setup directory not found:** check `PATH`,
  `PISA_SETUP_DIR`, or the corresponding CLI options. Both are required even when
  the assembly command reuses existing XML.
- **Configuration or PISA data-file errors:** check `pisa_cfg_tmp`, its placeholder
  for the selected workflow, and the `srs/`, `molref/`, and `pisastore/` paths.
- **Output file cannot be created:** create the JSON parent directory and check
  write permissions, including the container user's access for Docker runs.
- **Missing chain identifiers:** assign alphanumeric chain IDs to all polymer and
  ligand entities in the input structure; the validated parsers raise
  `MissingChainIDError` for unassigned IDs.

CLI argument errors are reported by argparse; PISA subprocess failures and
parsing/validation exceptions are propagated. The current application does not
define a dedicated exit status `9` for missing files.

## Versioning

We use [SemVer](https://semver.org) for versioning.

## Authors

- [Grisell Diaz Leines](https://github.com/grisell) - Lead developer
- [Stephen Anyango](https://github.com/otienoanyango) - Review and productionising
- [Mihaly Varadi](https://github.com/mvaradi) - Review and management

See all contributors [here](https://github.com/PDBe-KB/pisa-analysis/graphs/contributors).

## License

Apache License 2.0. See [LICENSE](LICENSE).
