import logging
import os
import subprocess
import tempfile
from typing import Tuple
from time import time

from pisa_utils.utils import create_pisa_config


def run_pisalite(
    input_cif, xml_output_dir, pisa_binary, pisa_setup_dir
) -> Tuple[str, str]:
    """Runs PISA Lite to determine interfaces, and returns XML files describing the assembly and interfaces.

    Expects PISA_

    Args:
        input_cif (str): Path to input CIF file.
        xml_output_dir (str): Path to output directory for XML files.
        pisa_binary (str, optional): Path to PISA binary. Defaults to "pisa".
        pisa_setup_dir (str, optional): Path to PISA setup directory.
            Checks PISA_SETUP_DIR environment variable if not provided. Defaults to None.

    Returns:
        Tuple[str, str]: Paths to XML files describing the assembly and interfaces.
    """
    
    
    with tempfile.TemporaryDirectory() as temp_dir:
        start = time()
        os.makedirs(xml_output_dir, exist_ok=True)
        cfg_file = create_pisa_config(temp_dir, pisa_setup_dir)
        
        session_name = "XXX"
        
        xml_interfaces_file = os.path.join(xml_output_dir, "interfaces.xml")
        xml_assembly_file = os.path.join(xml_output_dir, "assembly.xml")

        logging.info(f"Running {pisa_binary} on {input_cif}")
        subprocess.run(
            [pisa_binary, session_name, "-analyse", input_cif, cfg_file],
            check=True,
        )
        with open(xml_interfaces_file, "w") as f:
            subprocess.run(
                [pisa_binary, session_name, "-xml", "interfaces", cfg_file],
                stdout=f,
                check=True,
            )

        with open(xml_assembly_file, "w") as f:
            subprocess.run(
                [pisa_binary, session_name, "-xml", "assemblies", cfg_file],
                stdout=f,
                check=True,
            )
        if not (
            os.path.getsize(xml_interfaces_file) > 0
            and os.path.getsize(xml_assembly_file) > 0
        ):
            raise Exception(
                f"One or both XML files are empty: {xml_interfaces_file}, {xml_assembly_file}"
            )
        logging.info(f"XML files: {xml_assembly_file}, {xml_interfaces_file}")
        totaltime=time()-start
        
        logging.info(f"Finished analysis of interfaces in {totaltime} seconds")
    return xml_assembly_file, xml_interfaces_file
