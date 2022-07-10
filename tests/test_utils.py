import os
from unittest import TestCase

from pisa_utils.utils import parse_xml_file, read_uniprot_info


class TestUtils(TestCase):
    def test_xml_parser_valid(self):
        # Test that the helper function returns data when XML input exists
        self.assertTrue(parse_xml_file("./tests/data/basic_mock.xml"))

    def test_xml_parser_invalid(self):
        # Test that the helper function returns None when XML input does not exist
        self.assertIsNone(parse_xml_file("./tests/data/invalid.xml"))

    def test_read_uniprot_info_valid(self):
        """
        Test that the function works as intended
        """
        path_to_updated_cif = os.path.join(".", "tests", "data")
        result = read_uniprot_info("5", "2", "N", "ALA", "6nxr", path_to_updated_cif)
        self.assertEqual(result, (("P48491", "2")))
        result = read_uniprot_info("6", "3", "CA", "ARG", "6nxr", path_to_updated_cif)
        self.assertEqual(result, (("P48491", "3")))

    def test_read_uniprot_info_valid_with_negative_number(self):
        """
        Test that the function handles negative sequence number
        """
        path_to_updated_cif = os.path.join(".", "tests", "data")
        result = read_uniprot_info("-1", "3", "O", "ARG", "6nxr", path_to_updated_cif)
        self.assertEqual(result, ("P48491", "3"))

    def test_read_uniprot_info_no_atom(self):
        """
        Test that the function returns None when it can't find the atom
        """
        path_to_updated_cif = os.path.join(".", "tests", "data")
        result = read_uniprot_info("6", "3", "CA", "GLY", "6nxr", path_to_updated_cif)
        self.assertEqual(result, ((None, None)))

    def test_read_uniprot_info_no_mapping(self):
        """
        Test that the function returns None when it can't find UniProt mapping
        """
        path_to_updated_cif = os.path.join(".", "tests", "data")
        result = read_uniprot_info("6", "3", "CB", "ARG", "6nxr", path_to_updated_cif)
        self.assertEqual(result, ((None, None)))
