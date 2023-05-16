import os
import os.path
from unittest import TestCase
import tempfile

from gemmi import cif

from pisa_utils.utils import create_pisa_config, parse_xml_file, read_uniprot_info


class TestUtils(TestCase):
    def test_xml_parser_valid(self):
        # Test that the helper function returns data when XML input exists
        self.assertTrue(parse_xml_file("./tests/data/basic_mock.xml"))

    def test_xml_parser_invalid(self):
        # Test that the helper function returns None when XML input does not exist
        with self.assertRaises(FileNotFoundError):
            parse_xml_file("./tests/data/invalid.xml")

    def test_create_pisa_config(self):
        # test if function creates pisa configuration file
        dataroot = os.path.join(".", "tests", "data")
        setup_dir = os.path.join(".", "tests", "data")
        result = create_pisa_config(dataroot, setup_dir)
        pisa_cfg_exists = os.path.exists(result)

        self.assertTrue(pisa_cfg_exists)

    def test_read_uniprot_info_valid(self):
        """
        Test that the function works as intended
        """

        path_to_updated_cif = os.path.join(".", "tests", "data", "6nxr_updated.cif")
        doc = cif.read(path_to_updated_cif)
        updated_cif_block = doc.sole_block()

        result = read_uniprot_info("5", "2", "N", "ALA", updated_cif_block)
        self.assertEqual(result, (("P48491", "2")))
        result = read_uniprot_info("6", "3", "CA", "ARG", updated_cif_block)
        self.assertEqual(result, (("P48491", "3")))

    def test_read_uniprot_info_valid_with_negative_number(self):
        """
        Test that the function returns none with negative sequence number
        """
        path_to_updated_cif = os.path.join(".", "tests", "data", "6nxr_updated.cif")
        doc = cif.read(path_to_updated_cif)
        updated_cif_block = doc.sole_block()

        result = read_uniprot_info("-1", "3", "O", "ARG", updated_cif_block)
        self.assertEqual(result, (None, None))

    def test_read_uniprot_info_no_atom(self):
        """
        Test that the function returns None when it can't find the atom
        """
        path_to_updated_cif = os.path.join(".", "tests", "data", "6nxr_updated.cif")
        doc = cif.read(path_to_updated_cif)
        updated_cif_block = doc.sole_block()
        result = read_uniprot_info("6", "3", "CA", "GLY", updated_cif_block)
        self.assertEqual(result, ((None, None)))

    def test_read_uniprot_info_no_mapping(self):
        """
        Test that the function returns None when it can't find UniProt mapping
        """
        path_to_updated_cif = os.path.join(".", "tests", "data", "6nxr_updated.cif")
        doc = cif.read(path_to_updated_cif)
        updated_cif_block = doc.sole_block()
        result = read_uniprot_info("4", "3", "CB", "ARG", updated_cif_block)
        self.assertEqual(result, ((None, None)))
