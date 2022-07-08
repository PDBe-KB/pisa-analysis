import unittest
from unittest import TestCase
import os
import xmltodict

from pisa_utils import analyse
test_data_path= os.path.join('tests', 'data')

class TestXmlFiles(TestCase):
        
    def test_xml_pisa_output(self):

        ref_interfaces_data_file = os.path.join(test_data_path,'test_interfaces.xml')
        ref_assembly_data_file  = os.path.join(test_data_path,'test_assembly.xml')
        
        with open(ref_interfaces_data_file,'r') as test_interfaces_file:
            xml_interfaces=test_interfaces_file.read()
            
        ref_interfaces_dict = xmltodict.parse(xml_interfaces)
        
        with open(ref_assembly_data_file) as test_assembly_file:
            xml_assembly=test_assembly_file.read()
            
        ref_assembly_dict = xmltodict.parse(xml_assembly)

        pisa_config_file=os.path.join(test_data_path,'pisa.cfg')

        ap= analyse.AnalysePisa(pdbid_id="6nxr", assembly_id="1", pisa_config = pisa_config_file, output_dir=test_data_path, force=True, result_json_file="output.json", input_dir=test_data_path)
        
        ap.run_process()

        output_interfaces_data_file = os.path.join(test_data_path,'interfaces.xml')
        output_assembly_data_file = os.path.join(test_data_path,'assembly.xml')
        
        with open(output_interfaces_data_file,'r') as output_interfaces_file:
            out_interfaces_xml = output_interfaces_file.read()

        out_interfaces_dict=xmltodict.parse(out_interfaces_xml)
        
        with open(output_assembly_data_file,'r') as output_assembly_file:
            out_assembly_xml=output_assembly_file.read()

        out_assembly_dict=xmltodict.parse(out_assembly_xml)

        self.assertEqual(ref_interfaces_dict, out_interfaces_dict)
        self.assertEqual(ref_assembly_dict, out_assembly_dict)
                                                 
        
if __name__ == '__main__':
    unittest.main()

        

