import unittest
from unittest import TestCase
import os
import xmlschema

from pisa_to_json import main
test_data_path= os.path.join('tests', 'data')

class TestXmlFiles(TestCase):
        
    def test_validate_xml(self):

        schema_interfaces = os.path.join(test_data_path,'interfaces_schema.xsd')
        schema_assembly = os.path.join(test_data_path,'assembly_schema.xsd')

        path_pisa_config=os.path.join(test_data_path,'pisa.cfg')
        
        ap= main.AnalysePisa(pdbid_id="6nxr", assembly_id="1", pisa_config = path_pisa_config,output_dir=test_data_path,force=True,result_json_file="output.json",input_dir=test_data_path)
        
        ap.run_process()

        output_interfaces_xml = os.path.join(test_data_path,'interfaces.xml')
        output_assembly_xml = os.path.join(test_data_path,'assembly.xml')
        
        xmlschema.validate(output_interfaces_xml,schema_interfaces)
        xmlschema.validate(output_assembly_xml,schema_assembly)
                                                 
        
if __name__ == '__main__':
    unittest.main()

        

