import os
from unittest import TestCase

from pisa_utils.utils import parse_xml_file, read_uniprot_info


class TestUtils(TestCase):
    def test_xml_parser(self):
        # Test that the helper function returns data when XML input exists
        self.assertTrue(parse_xml_file("./tests/data/mocks/basic_mock.xml"))
        # Test that the helper function returns None when XML input does not exist
        self.assertIsNone(parse_xml_file("./tests/data/mocks/invalid.xml"))

    def test_read_uniprot_info(self):
        path_to_updated_cif = os.path.join(".", "tests", "data", "mocks")
        result = read_uniprot_info("5", "2", "N", "ALA", "6nxr", path_to_updated_cif)
        self.assertEqual(result, (("P48491", "2")))
        result = read_uniprot_info("6", "3", "CA", "ARG", "6nxr", path_to_updated_cif)
        self.assertEqual(result, (("P48491", "3")))
        result = read_uniprot_info("6", "3", "CA", "GLY", "6nxr", path_to_updated_cif)
        self.assertEqual(result, ((None, None)))
