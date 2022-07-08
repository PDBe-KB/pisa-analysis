import os
import unittest
from unittest import TestCase

import xmlschema

from pisa_utils import analyze

test_data_path = os.path.join("tests", "data")


class TestXmlFiles(TestCase):
    def test_validate_xml(self):

        schema_interfaces = os.path.join(test_data_path, "interfaces_schema.xsd")
        schema_assembly = os.path.join(test_data_path, "assembly_schema.xsd")

        ap = analyze.AnalysePisa(
            pdb_id="6nxr",
            assembly_id="1",
            output_dir=test_data_path,
            result_json_file="output.json",
            input_dir=test_data_path,
            input_updated_cif=None,
            input_cif_file=None,
        )

        ap.process_pisa_xml()
        ap.set_results()

        output_interfaces_xml = os.path.join(test_data_path, "interfaces.xml")
        output_assembly_xml = os.path.join(test_data_path, "assembly.xml")

        xmlschema.validate(output_interfaces_xml, schema_interfaces)
        xmlschema.validate(output_assembly_xml, schema_assembly)


if __name__ == "__main__":
    unittest.main()
