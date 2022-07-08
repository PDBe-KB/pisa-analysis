import logging
import xml.etree.ElementTree as ET


def parse_xml_file(xml_file):
    root = None
    try:
        logging.debug('parsing: {}'.format(xml_file))
        tree = ET.parse(xml_file)
        root = tree.getroot()
    except Exception as e:
        logging.error('invalid xml file: {}'.format(xml_file))
        logging.error(e)
    return root