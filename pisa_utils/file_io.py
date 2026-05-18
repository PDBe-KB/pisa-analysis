import json
import os
import shutil
import xml.etree.ElementTree as ET

from pisa_utils.field_handlers import logger


import gzip


def open_compressed(file_path: str, compressed: bool = False, mode: str = "rt"):
    """
    Open a file that may be compressed with gzip.

    :param file_path: Path to the file to open.
    :type file_path: str
    :param mode: Mode to open the file in. Defaults to "rt" (read text).
    :type mode: str, optional
    :return: File object.
    :rtype: file object
    """

    if compressed:
        if not file_path.endswith(".gz"):
            logger.error(f"Expected a .gz file for compressed input, got: {file_path}")
            raise ValueError(
                f"Expected a .gz file for compressed input, got: {file_path}"
            )

        logger.debug(f"Opening compressed file: {file_path}")
        return gzip.open(file_path, mode)

    else:
        logger.debug(f"Opening uncompressed file: {file_path}")
        return open(file_path, mode)


def save_json(data: dict, output_path: str, compressed: bool = False) -> None:
    """
    Save a dictionary to a JSON file, with optional gzip compression.

    :param data: Data to save to JSON.
    :type data: dict
    :param output_path: Path to the output JSON file.
    :type output_path: str
    :param compressed: Whether to compress the output file with gzip. Defaults to
        False.
    :type compressed: bool, optional
    """

    if compressed and not output_path.endswith(".gz"):
        output_path += ".gz"

    if compressed:
        with gzip.open(output_path, "wt") as json_file:
            json.dump(data, json_file, indent=4)
    else:
        with open(output_path, "w") as json_file:
            json.dump(data, json_file, indent=4)

    logger.info(f"JSON file written successfully: {output_path}")


def compress_existing_file(file_path: str) -> None:
    """
    Compress a file with gzip. Defaults to maximum compression level (-9) and
    overwrites the original file.

    :param file_path: Path to the file to compress.
    :type file_path: str
    :return: Path to the compressed file.
    :rtype: str
    """

    compressed_file_path = file_path + ".gz"

    with open(file_path, "rb") as f_in:
        with gzip.open(compressed_file_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    os.remove(file_path)
    logger.debug(f"Compressed file created at: {compressed_file_path}")


def create_pisa_config(dataroot, setup_dir, placeholder_dataroot="path_dataroot"):
    """
    Function creates a configuration file to run pisa

    :param dataroot : path to output pisa data and session
    :param setup_dir : path to setup directory of pisa
    :return : None

    """
    outputname = os.path.join(dataroot, "pisa.cfg")
    with open(os.path.join(setup_dir, "pisa_cfg_tmp")) as infile:
        with open(outputname, "w") as outfile:
            for line in infile:
                line = line.replace(placeholder_dataroot, dataroot).replace(
                    "path_to_setup", setup_dir
                )
                outfile.write(line)
    return outputname


def parse_xml_file(xml_file):
    """Parse XML file.

    Args:
        xml_file (str): Path to XML file.

    Returns:
        xml.etree.ElementTree.Element: Root element of XML file.

    """
    tree = ET.parse(xml_file)
    return tree.getroot()
