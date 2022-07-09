from pisa_utils.utils import read_uniprot_info


def get_molecules_dictionary(molecules):
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


def get_bond_dictionary(bonds, bondtype, pdb_id, input_updated_cif):
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

        uniprot_info_1 = read_uniprot_info(
            label_seqnum_1,
            seqnum_1,
            atname_1,
            res_1,
            pdb_id,
            input_updated_cif,
        )
        unp_acc_1 = uniprot_info_1[0]
        unp_num_1 = uniprot_info_1[1]

        uniprot_info_2 = read_uniprot_info(
            label_seqnum_2,
            seqnum_2,
            atname_2,
            res_2,
            pdb_id,
            input_updated_cif,
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
