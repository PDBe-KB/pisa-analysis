import logging
import os
from datetime import datetime
import subprocess as sub

logger = logging.getLogger()


def run_pisalite(session_name,
                 input_cif,
                 cfg_input,
                 output_dir,
                 pisa_binary):
    """
    Runs pisa-lite to analyse interfaces in assembly file and obtain xml files

    :param session_name: type str - used to identify sessions
    :param input_cif: type str - path to input cif file
    :param cfg_input: type str - path to configuration file
    :param output_dir: type str - path to output folder
    :param pisa_binary: type str - path to pisa binary
    :return: None
    """

    start = datetime.now()
    logging.info('starting Pisa on {}'.format(input_cif))

    binary = os.path.join(pisa_binary, "pisa") if pisa_binary else "pisa"

    xml_interfaces_file = open(os.path.join(output_dir, "interfaces.xml"), "w")
    xml_assembly_file = open(os.path.join(output_dir, "assembly.xml"), "w")

    sub.run([binary, session_name, "-analyse", input_cif, cfg_input])
    sub.run([binary, session_name, "-xml", "interfaces", cfg_input], stdout=xml_interfaces_file)
    sub.run([binary, session_name, "-xml", "assemblies", cfg_input], stdout=xml_assembly_file)

    xml_interfaces_file.close()
    xml_assembly_file.close()

    end = datetime.now()
    logging.info('finished Pisa analysis')
    print('finished Pisa analysis')
    time_taken = end - start
    time_taken_str = str(time_taken)
    logging.info('time taken {}'.format(time_taken_str))
    print('time taken {}'.format(time_taken_str))
    # TODO: RETURN TRUE IF THERE WAS REASONABLE TIME BETWEEN START AND END
    # TODO: OTHERWISE RETURN FALSE AND LOG A WARNING
