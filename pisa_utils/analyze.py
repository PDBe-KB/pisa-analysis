import json
import logging
import os

from gemmi import cif

from pisa_utils.utils import parse_xml_file

logger = logging.getLogger()


class AnalysePisa:
    """
    This class is responsible for parsing XML files generated by pisa, and
    converting their content to JSON format. It also add UniProt numbering
    using updated mmCIF file
    """

    def __init__(
        self,
        pdb_id,
        assembly_id,
        input_dir,
        input_updated_cif,
        input_cif_file,
        output_dir,
        result_json_file,
    ):
        self.pdb_id = pdb_id
        self.assembly_id = assembly_id
        self.input_dir = input_dir if input_dir else None
        self.input_updated_cif = input_updated_cif if input_updated_cif else output_dir
        self.input_cif_file = input_cif_file if input_cif_file else None
        self.output_dir = output_dir if output_dir else None
        self.result_json_file = result_json_file if result_json_file else None
        self.results = {}
        self.interfaces_results = None

    def _read_uniprot_info(self, int_lab_seqnum, int_seqnum, int_atname, int_resname):
        """
        Function reads atom's UniProt accession numbers and sequence numbers
        from updated cif file

        :param int_lab_seqnum: type str - atom label sequence id
        :param int_seqnum: type str - atom sequence number
        :param int_atname: type str - atom name
        :param int_resname: type str - atom residue name
        :return: (atom uniprot accession number, atom uniprot sequence number)
        """

        unp_acc = None
        unp_num = None

        # updated cif file path
        path = os.path.join(
            self.input_updated_cif, "{}_updated.cif".format(self.pdb_id)
        )

        # Reading UniProt acc and seq numbers in updated cif file:
        doc = cif.read_file(path)
        block = doc.sole_block()
        label_seq_id = block.find_loop("_atom_site.label_seq_id")
        atom_name = block.find_loop("_atom_site.label_atom_id")
        auth_seq_id = block.find_loop("_atom_site.auth_seq_id")
        db_name = block.find_loop("_atom_site.pdbx_sifts_xref_db_name")
        db_acc = block.find_loop("_atom_site.pdbx_sifts_xref_db_acc")
        db_num = block.find_loop("_atom_site.pdbx_sifts_xref_db_num")
        res_name = block.find_loop("_atom_site.label_comp_id")

        # counts if atom is not found in updated cif file
        n = 0

        # if sequence identifier read in pisa-lite is not available,
        # replace with a dot
        if int(int_lab_seqnum) < 0:
            int_lab_seqnum = "."

        # Search atom in updated cif file and read uniprot numbers
        for (labseqnum, seqnum, name, resname, dbname, dbacc, dbnum) in zip(
            label_seq_id, auth_seq_id, atom_name, res_name, db_name, db_acc, db_num
        ):
            if (
                labseqnum == int_lab_seqnum.strip()
                and seqnum == int_seqnum.strip()
                and name == int_atname.strip()
                and resname == int_resname.strip()
            ):
                n += 1
                if dbname == "UNP":
                    unp_acc = dbacc
                    unp_num = dbnum

                    return unp_acc, unp_num
                else:
                    logging.debug("No UNP numbers found for atom:")
                    logging.debug(
                        "name {},label_seq_id {},seq_num {}, residue {}".format(
                            int_atname, int_lab_seqnum, int_seqnum, int_resname
                        )
                    )
                    unp_acc = None
                    unp_num = None
                    return unp_acc, unp_num
        if n == 0:  # If atom was not found in updated cif file, return message
            logging.debug("atom not found in updated cif file:")
            logging.debug(
                "name {},label_seq_id {},seq_num {}, residue {}".format(
                    int_atname, int_lab_seqnum, int_seqnum, int_resname
                )
            )
            return unp_acc, unp_num

    def _get_bond_dictionary(self, bonds, bondtype):
        """
        Creates bond dictionary

        :param bonds: type dict - xml data dictionary
        :param bondtype: type str - bond interaction type
        :return: type dict - bonds dictionary
        """

        atom_site1_chains = []
        atom_site1_residues = []
        atom_site1_label_asym_ids = []
        atom_site1_orig_label_asym_ids = []
        atom_site1_unp_nums = []
        atom_site1_unp_accs = []
        atom_site1_seq_nums = []
        atom_site1_label_seq_ids = []
        atom_site1_label_atom_ids = []
        atom_site1_inscodes = []

        atom_site2_chains = []
        atom_site2_residues = []
        atom_site2_label_asym_ids = []
        atom_site2_orig_label_asym_ids = []
        atom_site2_unp_nums = []
        atom_site2_unp_accs = []
        atom_site2_seq_nums = []
        atom_site2_label_seq_ids = []
        atom_site2_label_atom_ids = []
        atom_site2_inscodes = []

        bond_types = []
        bonds_distances = []

        for bond in bonds:
            chain_1 = bond.find("chain-1").text
            label_asym_id_1 = bond.find("label_asym_id-1").text
            orig_label_asym_id_1 = bond.find("orig_label_asym_id-1").text
            res_1 = bond.find("res-1").text
            seqnum_1 = bond.find("seqnum-1").text
            label_seqnum_1 = bond.find("label_seqnum-1").text
            inscode_1 = bond.find("inscode-1").text
            atname_1 = bond.find("atname-1").text
            chain_2 = bond.find("chain-2").text
            label_asym_id_2 = bond.find("label_asym_id-2").text
            orig_label_asym_id_2 = bond.find("orig_label_asym_id-2").text
            res_2 = bond.find("res-2").text
            seqnum_2 = bond.find("seqnum-2").text
            label_seqnum_2 = bond.find("label_seqnum-2").text
            inscode_2 = bond.find("inscode-2").text
            atname_2 = bond.find("atname-2").text
            distance = bond.find("dist").text
            dist = round(float(distance), 2)
            bond_type = bondtype

            # *** Read uniprot accession and sequence numbers for
            # atoms in bonds, from updated cif file

            uniprot_info_1 = self._read_uniprot_info(
                label_seqnum_1, seqnum_1, atname_1, res_1
            )
            unp_acc_1 = uniprot_info_1[0]
            unp_num_1 = uniprot_info_1[1]

            uniprot_info_2 = self._read_uniprot_info(
                label_seqnum_2, seqnum_2, atname_2, res_2
            )
            unp_acc_2 = uniprot_info_2[0]
            unp_num_2 = uniprot_info_2[1]

            atom_site1_chains.append(chain_1)
            atom_site1_residues.append(res_1)
            atom_site1_label_asym_ids.append(label_asym_id_1)
            atom_site1_orig_label_asym_ids.append(orig_label_asym_id_1)
            atom_site1_unp_nums.append(unp_num_1)
            atom_site1_unp_accs.append(unp_acc_1)
            atom_site1_seq_nums.append(seqnum_1)
            atom_site1_label_seq_ids.append(label_seqnum_1)
            atom_site1_label_atom_ids.append(atname_1)
            atom_site1_inscodes.append(inscode_1)

            atom_site2_chains.append(chain_2)
            atom_site2_residues.append(res_2)
            atom_site2_label_asym_ids.append(label_asym_id_2)
            atom_site2_orig_label_asym_ids.append(orig_label_asym_id_2)
            atom_site2_unp_nums.append(unp_num_2)
            atom_site2_unp_accs.append(unp_acc_2)
            atom_site2_seq_nums.append(seqnum_2)
            atom_site2_label_seq_ids.append(label_seqnum_2)
            atom_site2_label_atom_ids.append(atname_2)
            atom_site2_inscodes.append(inscode_2)

            bonds_distances.append(dist)
            bond_types.append(bond_type)

        bond_dict = {
            "bond_distances": bonds_distances,
            "atom_site_1_chains": atom_site1_chains,
            "atom_site_1_residues": atom_site1_residues,
            "atom_site_1_label_asym_ids": atom_site1_label_asym_ids,
            "atom_site_1_orig_label_asym_ids": atom_site1_orig_label_asym_ids,
            "atom_site_1_unp_accs": atom_site1_unp_accs,
            "atom_site_1_unp_nums": atom_site1_unp_nums,
            "atom_site_1_seq_nums": atom_site1_seq_nums,
            "atom_site_1_label_seq_ids": atom_site1_label_seq_ids,
            "atom_site_1_label_atom_ids": atom_site1_label_atom_ids,
            "atom_site_1_inscodes": atom_site1_inscodes,
            "atom_site_2_chains": atom_site2_chains,
            "atom_site_2_residues": atom_site2_residues,
            "atom_site_2_label_asym_ids": atom_site2_label_asym_ids,
            "atom_site_2_orig_label_asym_ids": atom_site2_orig_label_asym_ids,
            "atom_site_2_unp_accs": atom_site2_unp_accs,
            "atom_site_2_unp_nums": atom_site2_unp_nums,
            "atom_site_2_seq_nums": atom_site2_seq_nums,
            "atom_site_2_label_seq_ids": atom_site2_label_seq_ids,
            "atom_site_2_label_atom_ids": atom_site2_label_atom_ids,
            "atom_site_2_inscodes": atom_site2_inscodes,
        }

        return bond_dict

    def _get_molecules_dictionary(self, molecules):
        """
        Function creates molecules and residues dictionaries from xml data

        :param molecules: xml data for interfaces molecues and residues
        :return: molecules_dicts: type dict - molecules dictionary
        :return: interface_residues_count : type int - residue count
        :return: is_ligand: type bool - is molecule class a 'Ligand'?
        """

        is_ligand = False
        molecules_dicts = []

        for molecule in molecules:
            interface_residues_count = 0
            molecule_id = molecule.find("id").text
            molecule_class = molecule.find("class").text
            chain_id = molecule.find("chain_id").text
            residues_dicts = []

            if molecule_class in ["Ligand"]:
                is_ligand = True
            interface_residues = molecule.findall("residues/residue")
            residue_label_ids = []
            residue_sequence_numbers = []
            residue_label_sequence_numbers = []
            residue_ins_codes = []
            accessible_surface_areas = []
            buried_surface_areas = []
            solvation_energy_effects = []
            residue_bonds = []

            # Creating residues dictionaries
            for residue in interface_residues:
                residue_sernum = residue.find("ser_no").text
                residue_name = residue.find("name").text
                residue_seqnum = residue.find("seq_num").text
                residue_label_seq_num = residue.find("label_seq_num").text
                residue_asa = round(float(residue.find("asa").text), 2)
                residue_bsa = round(float(residue.find("bsa").text), 2)
                residue_solv_en = round(float(residue.find("solv_en").text), 2)
                residue_ins_code = residue.find("ins_code").text
                residue_bond = residue.find("bonds").text

                # Writing interface residues dictionary
                residue_dict = {
                    "residue_sernum": residue_sernum,
                    "residue_name": residue_name,
                    "residue_seqnum": residue_seqnum,
                    "residue_label_seq_num": residue_label_seq_num,
                    "residue_asa": residue_asa,
                    "residue_bsa": residue_bsa,
                    "residue_ins_code": residue_ins_code,
                    "residue_solv_en": residue_solv_en,
                    "residue_bond": residue_bond,
                }
                residues_dicts.append(residue_dict)
                residue_label_ids.append(residue_name)
                residue_sequence_numbers.append(residue_seqnum)
                residue_label_sequence_numbers.append(residue_label_seq_num)
                accessible_surface_areas.append(residue_asa)
                residue_ins_codes.append(residue_ins_code)
                residue_bonds.append(residue_bond)
                buried_surface_areas.append(residue_bsa)
                solvation_energy_effects.append(residue_solv_en)
                interface_residues_count += 1

            # Writing molecules dictionary
            molecule_dict = {
                "molecule_id": molecule_id,
                "molecule_class": molecule_class,
                "chain_id": chain_id,
                "residue_label_comp_ids": residue_label_ids,
                "residue_seq_ids": residue_sequence_numbers,
                "residue_label_seq_ids": residue_label_sequence_numbers,
                "residue_ins_codes": residue_ins_codes,
                "residue_bonds": residue_bonds,
                "solvation_energies": solvation_energy_effects,
                "accessible_surface_areas": accessible_surface_areas,
                "buried_surface_areas": buried_surface_areas,
            }
            molecules_dicts.append(molecule_dict)

            # if there is only one inteface residues,
            # discard interface as valid interface
            if len(interface_residues) == 1:
                is_ligand = True

        return molecules_dicts, interface_residues_count, is_ligand

    def process_pisa_xml(self):
        """
        Function writes assembly interfaces dictionaries
        :return: type dict - interfaces dictionary
        """

        interfaces_xml_file = os.path.join(self.output_dir, "interfaces.xml")
        assembly_xml_file = os.path.join(self.output_dir, "assembly.xml")
        result = {}

        if os.path.exists(interfaces_xml_file) and os.path.exists(assembly_xml_file):
            asroot = parse_xml_file(xml_file=assembly_xml_file)
            root = parse_xml_file(xml_file=interfaces_xml_file)
            if root and asroot:
                assembly_status = asroot.find("status").text
                assemblies = asroot.findall("asu_complex")

                # Assembly information
                for assem in assemblies:
                    assem_mmsize = assem.find("assembly/mmsize").text
                    assem_diss_energy = assem.find("assembly/diss_energy").text
                    assem_asa = assem.find("assembly/asa").text
                    assem_bsa = assem.find("assembly/bsa").text
                    assem_entropy = assem.find("assembly/entropy").text
                    assem_diss_area = assem.find("assembly/diss_area").text
                    assem_int_energy = assem.find("assembly/int_energy").text
                    assem_formula = assem.find("assembly/formula").text
                    assem_composition = assem.find("assembly/composition").text

                result["assembly_status"] = assembly_status

                # Round to two decimals some assembly properties
                assembly_mmsize = assem_mmsize
                assembly_diss_energy = round(float(assem_diss_energy), 2)
                assembly_asa = round(float(assem_asa), 2)
                assembly_bsa = round(float(assem_bsa), 2)
                assembly_entropy = round(float(assem_entropy), 2)
                assembly_diss_area = round(float(assem_diss_area), 2)
                assembly_int_energy = round(float(assem_int_energy), 2)
                assembly_formula = assem_formula
                assembly_composition = assem_composition

                # Create interfaces dictionaries
                status = root.find("status").text
                num_interfaces = root.find("n_interfaces").text
                interfaces = root.findall("interface")
                logging.debug("number of interfaces: {}".format(len(interfaces)))

                result["status"] = status
                result["num_interfaces"] = num_interfaces

                non_ligand_interface_count = 0

                for interface in interfaces:

                    # Interface General information
                    interface_id = interface.find("id").text
                    interface_area = round(float(interface.find("int_area").text), 2)
                    interface_solvation_energy = round(
                        float(interface.find("int_solv_en").text), 2
                    )
                    interface_stabilization_energy = round(
                        float(interface.find("stab_en").text), 2
                    )
                    p_value = round(float(interface.find("pvalue").text), 3)

                    # No. of bonds counted

                    n_h_bonds = int(interface.find("h-bonds/n_bonds").text)
                    n_ss_bonds = int(interface.find("ss-bonds/n_bonds").text)
                    n_covalent_bonds = int(interface.find("cov-bonds/n_bonds").text)
                    n_salt_bridges = int(interface.find("salt-bridges/n_bonds").text)
                    other_contacts = int(interface.find("other-bonds/n_bonds").text)

                    # Reading bonds
                    hbonds = interface.findall("h-bonds/bond")
                    sbridges = interface.findall("salt-bridges/bond")
                    covbonds = interface.findall("cov-bonds/bond")
                    ssbonds = interface.findall("ss-bonds/bond")
                    othbonds = interface.findall("other-bonds/bond")

                    # Writing bonds dictionaries
                    hbond_dict = self._get_bond_dictionary(hbonds, "H-bond")
                    sbridge_dict = self._get_bond_dictionary(sbridges, "salt-bridges")
                    covbond_dict = self._get_bond_dictionary(covbonds, "cov-bonds")
                    ssbond_dict = self._get_bond_dictionary(ssbonds, "ss-bonds")
                    othbond_dict = self._get_bond_dictionary(othbonds, "other-bond")
                    molecules = interface.findall("molecule")

                    (
                        molecules_dicts,
                        interface_residues_count,
                        is_ligand,
                    ) = self._get_molecules_dictionary(molecules)

                    if not is_ligand:
                        non_ligand_interface_count += 1

                        interface_dict = {
                            "interface_id": interface_id,
                            "interface_area": interface_area,
                            "solvation_energy": interface_solvation_energy,
                            "stabilization_energy": interface_stabilization_energy,
                            "p_value": p_value,
                            "number_interface_residues": interface_residues_count,
                            "number_hydrogen_bonds": n_h_bonds,
                            "number_covalent_bonds": n_covalent_bonds,
                            "number_disulfide_bonds": n_ss_bonds,
                            "number_salt_bridges": n_salt_bridges,
                            "number_other_bonds": other_contacts,
                            "hydrogen_bonds": hbond_dict,
                            "salt_bridges": sbridge_dict,
                            "disulfide_bonds": ssbond_dict,
                            "covalent_bonds": covbond_dict,
                            "other_bonds": othbond_dict,
                            "molecules": molecules_dicts,
                        }

                        # Append all dictionaries in 'result'

                        result.setdefault("id", []).append(interface_id)
                        result.setdefault("int_area", []).append(interface_area)
                        result.setdefault("interface_dicts", []).append(interface_dict)

                # Assembly information added to dictionary

                result["non_ligand_interface_count"] = non_ligand_interface_count
                result["assembly_mmsize"] = assembly_mmsize
                result["assembly_diss_energy"] = assembly_diss_energy
                result["assembly_asa"] = assembly_asa
                result["assembly_bsa"] = assembly_bsa
                result["assembly_entry"] = assembly_entropy
                result["assembly_diss_area"] = assembly_diss_area
                result["assembly_int_energy"] = assembly_int_energy
                result["assembly_formula"] = assembly_formula
                result["assembly_composition"] = assembly_composition

        return result

    def set_results(self):
        """
        Writes assembly dictionary

        Args:
        interfaces_results: interfaces dictionaries
        entry_type: type sting - entry type 'assembly'
        pdb_id : pdb entry
        assembly_id: assembly code

        Returns: None

        """
        if self.interfaces_results:
            overall = len(self.interfaces_results.get("id", []))

            interface_dicts = self.interfaces_results.get("interface_dicts", [])

            assem_dict = {
                "mmsize": self.interfaces_results.get("assembly_mmsize"),
                "dissociation_energy": self.interfaces_results.get(
                    "assembly_diss_energy"
                ),
                "accessible_surface_area": self.interfaces_results.get("assembly_asa"),
                "buried_surface_area": self.interfaces_results.get("assembly_bsa"),
                "entropy": self.interfaces_results.get("assembly_entropy"),
                "dissociation_area": self.interfaces_results.get("assembly_diss_area"),
                "solvation_energy_gain": self.interfaces_results.get(
                    "assembly_int_energy"
                ),
                "formula": self.interfaces_results.get("assembly_formula"),
                "composition": self.interfaces_results.get("assembly_composition"),
                "interface_count": overall,
                "interfaces": interface_dicts,
            }

            assembly_dictionary = {
                "pdb_id": self.pdb_id,
                "assembly_id": self.assembly_id,
                "pisa_version": "2.0",
                "assembly": assem_dict,
            }
            self.results.setdefault("PISA", assembly_dictionary)

    def save_to_json(self):
        """
        Dump the data into a JSON file
        :return:
        """
        # TODO - ADD TEST
        # TODO - REFACTOR
        # TODO - TRUE IF SUCCESS FALSE IF FAIL
        # TODO - ADD LOGGING

        if self.results:
            output_file = (
                os.path.join(self.output_dir, self.result_json_file)
                if self.result_json_file
                else os.path.join(
                    self.output_dir,
                    "{}-assembly{}.json".format(self.pdb_id, self.assembly_id),
                )
            )
            with open(output_file, "w") as out_file:
                json.dump(self.results, out_file)
