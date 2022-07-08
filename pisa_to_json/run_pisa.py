import logging
import os
import argparse
import sys
from datetime import datetime
import subprocess as sub

logger = logging.getLogger()

class RunPisa:
    
    def run_pisalite(self, session_name,input_cif, cfg_input,output_dir, pisa_binary):

        #************ Runs pisa-lite to analyse interfaces in assembly file and obtain xml files***********
        
        start = datetime.now()
        logging.info('starting Pisa on {}'.format(input_cif))
        
        binary=os.path.join(pisa_binary,"pisa") if pisa_binary else "pisa"
        xml_interfaces_file=open(os.path.join(output_dir,"interfaces.xml"),"w")
        xml_assembly_file=open(os.path.join(output_dir,"assembly.xml"),"w")
        proc1=sub.run([binary,session_name,"-analyse",input_cif,cfg_input])
        proc2=sub.run([binary,session_name,"-xml","interfaces",cfg_input],stdout=xml_interfaces_file)
        proc3=sub.run([binary,session_name,"-xml","assemblies",cfg_input],stdout=xml_assembly_file)
        #proc4=sub.run(["mv","session_name","-xml","assemblies",cfg_input])

        xml_interfaces_file.close()
        xml_assembly_file.close()

        end = datetime.now()
        logging.info('finished Pisa analysis')
        print('finished Pisa analysis')
        time_taken = end-start
        time_taken_str = str(time_taken)
        logging.info('time taken {}'.format(time_taken_str))
        print('time taken {}'.format(time_taken_str))
        
        return None
            
    def run_process(self, session_name,input_cif,cfg_input,output_dir,pisa_binary):

        #******** Running pisa-lite process (stands alone) on assembly structure ****
        
        logging.debug('session directory: {}'.format(output_dir))
        #binary=os.path.join(pisa_binary,"pisa") if pisa_binary else "pisa"
        return self.run_pisalite(session_name,input_cif=input_cif,cfg_input=cfg_input,output_dir=output_dir,pisa_binary=pisa_binary)
        


if '__main__' in __name__:

    #************** Running analysis of assemblies with pisa-lite and obtain xml files***************
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', help='debugging', action='store_const', dest='loglevel', const=logging.DEBUG,
                        default=logging.INFO)
    parser.add_argument('--pdb_id', help='pdb id and session name', type=str, required=True)
    parser.add_argument('--input_cif', help='input mmcif file', type=str, required=True)
    parser.add_argument('--cfg_input', help='pisa configuration file', type=str, required=True)
    parser.add_argument('--output_dir', help='output pisa interface and assembly xml files', type=str, required=True)
    parser.add_argument('--pisa_binary', help='pisa configuration file', type=str)

    args = parser.parse_args()
    logger.setLevel(args.loglevel)

    
    rp = RunPisa()
    ok = rp.run_process(session_name=args.pdb_id,input_cif=args.input_cif,cfg_input=args.cfg_input,output_dir=args.output_dir,pisa_binary=args.pisa_binary)
    
    #ok = rp.run_process(input_cif=args.input_cif, interaction_xml=args.output_xml)
#   if not ok:
#        logging.error('failed pisa analysis')
