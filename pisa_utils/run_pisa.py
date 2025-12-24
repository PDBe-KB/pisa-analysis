import logging
import os
import subprocess
import tempfile
from time import time
from typing import Optional, Tuple

from pisa_utils.constants import SUBDIR_EXTENDED_DATA
from pisa_utils.models.models import DataModes, LigandProcessingMode
from pisa_utils.utils import create_pisa_config


def run_pisalite(
    input_cif, xml_output_dir, pisa_binary, pisa_setup_dir
) -> Tuple[str, str]:
    """Runs PISA to determine interfaces, and returns XML files describing the
    assembly and interfaces.

    Expects PISA_

    Args:
        input_cif (str): Path to input CIF file.
        xml_output_dir (str): Path to output directory for XML files.
        pisa_binary (str, optional): Path to PISA binary. Defaults to "pisa".
        pisa_setup_dir (str, optional): Path to PISA setup directory.
            Checks PISA_SETUP_DIR environment variable if not provided.
            Defaults to None.

    Returns:
        Tuple[str, str]: Paths to XML files describing the assembly and interfaces.
    """

    with tempfile.TemporaryDirectory() as temp_dir:

        start = time()
        cfg_file = create_pisa_config(temp_dir, pisa_setup_dir)

        xml_assembly_file, xml_interfaces_file = _run_pisa_command(
            input_cif=input_cif,
            xml_output_dir=xml_output_dir,
            cfg_file=cfg_file,
            pisa_binary=pisa_binary,
        )

        totaltime = time() - start

        logging.info(f"Finished analysis of interfaces in {totaltime} seconds")

    return xml_assembly_file, xml_interfaces_file


def run_pisa_service(
    input_cif: str,
    xml_output_dir: str,
    pisa_binary: str,
    pisa_setup_dir: str,
    asis: bool = False,
    lig_proc_mode: LigandProcessingMode = None,
    exclude_ligands: list[str] = None,
    extended_data: bool = False,
) -> Tuple[str, str, Optional[str], Optional[str], Optional[str]]:
    """
    Wrapper to run PISA with most available options. Intended for internal use at the
    PDBe.

    :param input_cif: Path to input CIF file. Full crystallographic mmCIF is intended,
        but other mmCIF files could be used.
    :type input_cif: str
    :param xml_output_dir: Path to output directory for XML files.
    :type xml_output_dir: str
    :param pisa_binary: Path to PISA binary.
    :type pisa_binary: str
    :param pisa_setup_dir: Path to PISA setup directory.
    :type pisa_setup_dir: str
    :param asis: When true, calculate interfaces only and PISA will NOT predict
        assemblies, defaults to False
    :type asis: bool, optional
    :param lig_proc_mode: Mode in which PISA will handle ligands. Options are "auto",
        "fixed" and "free, defaults to None
    :type lig_proc_mode: LigandProcessingMode, optional
    :param exclude_ligands: List of ligand CCDs to exclude from interface/assembly
        analysis, defaults to None
    :type exclude_ligands: list[str], optional
    :param extended_data: Extract extract information via PISA using the -list
        operation, defaults to False
    :type extended_data: bool, optional
    :return: Paths to XML files describing assemblies and interfaces, and optionally
        extended data files for assemblies, interfaces and monomers
    :rtype: Tuple[str, str, Optional[str], Optional[str], Optional[str]]
    """

    with tempfile.TemporaryDirectory() as temp_dir:

        cfg_file = create_pisa_config(temp_dir, pisa_setup_dir, "/data")

        session_name = input_cif.split("/")[-1].split(".")[0].replace("_", "")

        output_paths = _run_pisa_command(
            input_cif=input_cif,
            xml_output_dir=xml_output_dir,
            session_name=session_name,
            cfg_file=cfg_file,
            pisa_binary=pisa_binary,
            asis=asis,
            lig_proc_mode=lig_proc_mode,
            exclude_ligands=exclude_ligands,
            extended_data=extended_data,
        )

        logging.info("Finished -analysis, -xml and -list steps of PISA")

    return output_paths


def _run_pisa_command(
    input_cif: str,
    xml_output_dir: str,
    session_name: str = "XXX",
    cfg_file: str = "pisa.cfg",
    pisa_binary: str = "pisa",
    asis: bool = False,
    lig_proc_mode: LigandProcessingMode = None,
    exclude_ligands: list[str] = None,
    extended_data: bool = False,
) -> tuple[str, str, Optional[str], Optional[str], Optional[str]]:
    """Wrapper to runa PISA command with various options.

    :param input_cif: mmCIF (or PDB) file to analyse
    :type input_cif: str
    :param xml_output_dir: Directory to write XML output files
    :type xml_output_dir: str
    :param session_name: Placeholder name for the -analyse session, defaults to "XXX"
    :type session_name: str, optional
    :param cfg_file: Path to PISA configuration file, defaults to "pisa.cfg"
    :type cfg_file: str, optional
    :param pisa_binary: Path to pre-installed PISA binary, defaults to "pisa"
    :type pisa_binary: str, optional
    :param asis: When true, calculate interfaces only and PISA will NOT predict
        assemblies, defaults to False
    :type asis: bool, optional
    :param lig_proc_mode: Mode in which PISA will handle ligands. Options are "auto",
        "fixed" and "free, defaults to None
    :type lig_proc_mode: LigandProcessingMode, optional
    :param exclude_ligands: List of ligand CCDs to exclude from interface/assembly
        analysis, defaults to None
    :type exclude_ligands: list[str], optional
    :param extended_data: Extract extract information via PISA using the -list
        operation, defaults to False
    :type extended_data: bool, optional
    :raises Exception: If one or both XML output files are empty
    :return: Paths to XML files describing assemblies and interfaces, and optionally
        extended data files for assemblies, interfaces and monomers
    :rtype: tuple[str, str, Optional[str], Optional[str], Optional[str]]
    """

    os.makedirs(xml_output_dir, exist_ok=True)

    # Base command
    command = [
        pisa_binary,
        session_name,
        "-analyse",
        input_cif,
    ]

    # Run in as-is mode
    if asis:
        logging.info("Running analysis in as-is mode")
        command.append("--as-is")

    # Exclude ligands
    if exclude_ligands:
        logging.info(f"Excluding ligands from analysis: {exclude_ligands}")
        with open(f"{xml_output_dir}/agents.dat", "w") as excl_lig_file:
            for ligand in exclude_ligands:
                excl_lig_file.write(f'{ligand}, "tmp", formula\n')

        command.append("--lig-exclude='(agents)," + ",".join(exclude_ligands) + "' ")

    # Ligand processing mode
    if lig_proc_mode:
        logging.info(f"Running analysis with ligand processing mode: {lig_proc_mode}")
        command.append(f"--lig={lig_proc_mode}")

    # Config
    command.append(cfg_file)

    logging.info(f"Running PISA command: {' '.join(command)}")

    # Run PISA -analyse
    subprocess.run(
        command,
        check=True,
    )

    # Run PISA -xml
    xml_interfaces_file = os.path.join(xml_output_dir, "interfaces.xml")
    xml_assembly_file = os.path.join(xml_output_dir, "assembly.xml")

    with open(xml_interfaces_file, "w") as f:
        subprocess.run(
            [pisa_binary, session_name, "-xml", DataModes.INTERFACES, cfg_file],
            stdout=f,
            check=True,
        )
    with open(xml_assembly_file, "w") as f:
        subprocess.run(
            [pisa_binary, session_name, "-xml", DataModes.ASSEMBLIES, cfg_file],
            stdout=f,
            check=True,
        )

    if not (
        os.path.getsize(xml_interfaces_file) > 0
        and os.path.getsize(xml_assembly_file) > 0
    ):
        raise Exception(
            f"One or both XML files are empty: {xml_interfaces_file}, "
            f"{xml_assembly_file}"
        )

    logging.info(f"XML files: {xml_assembly_file}, {xml_interfaces_file}")
    output_data_paths = [xml_assembly_file, xml_interfaces_file]

    # Run PISA -list
    if extended_data:
        for data_type in (
            DataModes.ASSEMBLIES,
            DataModes.INTERFACES,
            DataModes.MONOMERS,
        ):
            logging.info(f"Generating extended data for {data_type}")
            path = os.path.join(xml_output_dir, SUBDIR_EXTENDED_DATA)
            os.makedirs(path, exist_ok=True)
            path = os.path.join(path, f"{data_type}_extended.txt")

            with open(path, "w") as f:
                subprocess.run(
                    [pisa_binary, session_name, "-list", data_type, cfg_file],
                    stdout=f,
                    check=True,
                )

            logging.info(f"Extended data file generated: {path}")
            output_data_paths.append(path)

    return tuple(output_data_paths)
