import os
import os.path
from unittest import TestCase
import tempfile

from pisa_utils.utils import parse_xml_file, read_uniprot_info, create_pisa_config
from gemmi import cif

class TestUtils(TestCase):
    def test_xml_parser_valid(self):
        # Test that the helper function returns data when XML input exists
        self.assertTrue(parse_xml_file("./tests/data/basic_mock.xml"))

    def test_xml_parser_invalid(self):
        # Test that the helper function returns None when XML input does not exist
        self.assertIsNone(parse_xml_file("./tests/data/invalid.xml"))


    def test_create_pisa_config(self):
        # test if function creates pisa configuration file 
        #dataroot = os.path.join(".", "tests", "data")
        setup_dir = os.path.join(".", "tests", "data")
        with tempfile.TemporaryDirectory() as tempdir:
            result = create_pisa_config(tempdir, setup_dir)
            pisa_cfg_exists= os.path.exists(result)
        
            self.assertTrue(pisa_cfg_exists)
        
    def test_read_uniprot_info_valid(self):
        """
        Test that the function works as intended
        """

        path_to_updated_cif = os.path.join(".", "tests", "data","6nxr_updated.cif")
        doc = cif.read(path_to_updated_cif)
        updated_cif_block = doc.sole_block()
        
        result = read_uniprot_info("5", "ALA", "6nxr", updated_cif_block)
        self.assertEqual(result, (("P48491", "2")))
        result = read_uniprot_info("6", "ARG", "6nxr",updated_cif_block)
        self.assertEqual(result, (("P48491", "3")))
    
    def test_read_uniprot_info_valid_with_negative_number(self):
        """
        Test that the function returns none with negative sequence number
        """
        path_to_updated_cif = os.path.join(".", "tests", "data","6nxr_updated.cif")
        doc = cif.read(path_to_updated_cif)
        updated_cif_block = doc.sole_block()

        result = read_uniprot_info("-1", "ARG", "6nxr", updated_cif_block)
        self.assertEqual(result, (None, None))
    
    def test_read_uniprot_info_no_residue(self):
        """
        Test that the function returns None when it can't find the atom
        """
        path_to_updated_cif = os.path.join(".", "tests", "data","6nxr_updated.cif")
        doc = cif.read(path_to_updated_cif)
        updated_cif_block = doc.sole_block()
        result = read_uniprot_info("6", "GLY", "6nxr", updated_cif_block)
        self.assertEqual(result, ((None, None)))

    def test_read_uniprot_info_no_mapping(self):
        """
        Test that the function returns None when it can't find UniProt mapping
        """
        path_to_updated_cif = os.path.join(".", "tests", "data","6nxr_updated.cif")
        doc = cif.read(path_to_updated_cif)
        updated_cif_block = doc.sole_block()
        result = read_uniprot_info("4", "ARG", "6nxr", updated_cif_block)
        self.assertEqual(result, ((None, None)))
