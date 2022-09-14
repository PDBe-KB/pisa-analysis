import logging
import os
import shutil
import subprocess as sub
import tempfile
from datetime import datetime

from pisa_utils.utils import create_pisa_config

logger = logging.getLogger()


def run_pisalite(input_cif, output_xml, pisa_binary):
    """
    Runs pisa-lite to analyse interfaces in assembly file and obtain xml files

    :param input_cif: type str - path to input cif file
    :param output_dir: type str - path to output folder

    :return: None
    """

    # Create pisa configuration file:

    pisa_set_up = os.environ["PISA_SETUP_DIR"]
    create_pisa_config(output_xml, pisa_set_up)

    cfg_input = os.path.join(output_xml, "pisa.cfg")
    pisa_output = os.path.join(output_xml, "pisa_XXX")

    start = datetime.now()
    logging.info("starting Pisa on {}".format(input_cif))

    # binary and session name

    binary = os.path.join(pisa_binary, "pisa") if pisa_binary else "pisa"
    print(binary)
    session_name = "XXX"

    xml_interfaces_file = open(os.path.join(output_xml, "interfaces.xml"), "w")
    xml_assembly_file = open(os.path.join(output_xml, "assembly.xml"), "w")

    sub.run([binary, session_name, "-analyse", input_cif, cfg_input])
    sub.run(
        [binary, session_name, "-xml", "interfaces", cfg_input],
        stdout=xml_interfaces_file,
    )
    sub.run(
        [binary, session_name, "-xml", "assemblies", cfg_input],
        stdout=xml_assembly_file,
    )

    xml_interfaces_file.close()
    xml_assembly_file.close()

    os.remove(cfg_input)
    shutil.rmtree(pisa_output)

    end = datetime.now()
    logging.info("finished Pisa analysis")

    time_taken = end - start
    time_taken_str = str(time_taken)
    logging.info("time taken {}".format(time_taken_str))

    # TODO: RETURN TRUE IF THERE WAS REASONABLE TIME BETWEEN START AND END
    # TODO: OTHERWISE RETURN FALSE AND LOG A WARNING
