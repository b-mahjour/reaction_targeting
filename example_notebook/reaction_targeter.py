#!/usr/bin/env python

import sys
import time
import json
import pandas as pd
import rdkit 
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')


def check_input_target_molecule(target_molecule_smiles):
    # checks to see if input target molecule smiles is valid
    
    mol = Chem.MolFromSmiles(target_molecule_smiles)
    if not mol:
        return None, "invalid target molecule smiles"
    try:
        Chem.SanitizeMol(mol)
    except:
        return None, "failed to sanitize smiles"
    return mol, ""


def check_input_target_substructure(target_substructure_smarts):
    # checks to see if input target molecule smiles is valid
    
    mol = Chem.MolFromSmarts(target_substructure_smarts)
    if not mol:
        return None, "invalid target substructure smarts"

    for atm in mol.GetAtoms():
        if atm.GetAtomMapNum() == 0:
            return None, "target substructure missing atom number label"
    
    return mol, ""


def check_number_of_substructure_hits(target_mol, substructure_mol):
    hits = target_mol.GetSubstructMatches(substructure_mol)
    
    if len(hits) == 0:
        return None, "no substructure hits in target"
    else:
        return hits, ""
    
    
def copy_and_map_atoms_to_target(target_mol, hit_atom_numbers, substructure_atom_map):
    copy = Chem.Mol(target_mol)
    
    for i, atm_num in enumerate(hit_atom_numbers):
        atm = copy.GetAtomWithIdx(atm_num)
        atm.SetAtomMapNum(substructure_atom_map[i])
    
    return copy, Chem.MolToSmiles(copy)


def _copy_atom(atom):
    new_atom = Chem.Atom(atom.GetSymbol())
    new_atom.SetFormalCharge(atom.GetFormalCharge())
    new_atom.SetAtomMapNum(atom.GetAtomMapNum())
    if atom.GetIsAromatic() and atom.GetSymbol() == 'N':
        new_atom.SetNumExplicitHs(atom.GetTotalNumHs())
    return new_atom


def disconnect_synthon(mapped_sm, dis_labels):
    mol = Chem.MolFromSmiles(mapped_sm)
    Chem.Kekulize(mol, clearAromaticFlags=False)
    Chem.SanitizeMol(mol)
    
    new_mol = Chem.RWMol(Chem.MolFromSmiles(''))
    for atom in mol.GetAtoms():
        new_atom = _copy_atom(atom)
        new_mol.AddAtom(new_atom)

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        a1idx = bond.GetBeginAtom().GetIdx()
        a2idx = bond.GetEndAtom().GetIdx()

        if [a1.GetAtomMapNum(), a2.GetAtomMapNum()] == dis_labels or \
            [a2.GetAtomMapNum(), a1.GetAtomMapNum()] == dis_labels:
            continue
        bt = bond.GetBondType()
        new_mol.AddBond(a1idx, a2idx, bt)

    # Chem.SanitizeMol(new_mol)
    transformed_smiles = Chem.MolToSmiles(new_mol).split(".")
    for tsm in transformed_smiles:
        if f":{str(dis_labels[0])}]" in tsm:
            synthon_a = tsm
        if f":{str(dis_labels[1])}]" in tsm:
            synthon_b = tsm
    synthon_a = synthon_a.replace("[n:1]", "[nH:1]")
    synthon_b = synthon_b.replace("[n:2]", "[nH:2]")
    return synthon_a, synthon_b


def remove_atom_labels(input_synthon_smiles):
    synthon_copy = Chem.MolFromSmiles(input_synthon_smiles)
    for i in synthon_copy.GetAtoms():
        i.SetAtomMapNum(0)
    return Chem.MolToSmiles(synthon_copy)


def check_enumerated_synthon_for_commercial_availability(cleansed_enumerated_synthon, commercial_dataset):
    for i,k in commercial_dataset.iterrows():
        if k["sanitized_smiles"] == cleansed_enumerated_synthon:
            return k.to_dict()
    return None


def add_fg_to_synthon(synthon_smiles, connection_label, building_block, building_block_label):
    mol = Chem.MolFromSmiles(synthon_smiles + "." + building_block)
    fg_obj = Chem.MolFromSmarts(building_block)
    
    new_mol = Chem.RWMol(Chem.MolFromSmiles(''))
    for atom in mol.GetAtoms():
        new_atom = _copy_atom(atom)
        new_mol.AddAtom(new_atom)

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        a1idx = bond.GetBeginAtom().GetIdx()
        a2idx = bond.GetEndAtom().GetIdx()

        bt = bond.GetBondType()
        new_mol.AddBond(a1idx, a2idx, bt)

    for atm in new_mol.GetAtoms():
        if atm.GetAtomMapNum() == connection_label:
            atom_to_add1 = atm.GetIdx()
        if atm.GetAtomMapNum() == building_block_label:
            atom_to_add2 = atm.GetIdx()

    new_mol.AddBond(atom_to_add1, atom_to_add2, rdkit.Chem.rdchem.BondType.SINGLE)
    enumerated_synthon = Chem.MolToSmiles(new_mol)
    if building_block_label != "hydrogen":
        enumerated_synthon = enumerated_synthon.replace("[nH:2]", "[n:2]")
    return enumerated_synthon


def check_synthon_for_commercial_analogs(synthon_smiles, connection_label, commercial_dataset, \
                                         building_blocks, building_block_label):
    enumerated_building_blocks = []
    enumerated_building_blocks.append(synthon_smiles)
    for i in building_blocks:
        enumerated_building_blocks.append(add_fg_to_synthon(synthon_smiles, connection_label, i, building_block_label))
    
    output_hits = []
    for i in enumerated_building_blocks:
        clean_smiles = remove_atom_labels(i)
        output_hits.append(check_enumerated_synthon_for_commercial_availability(clean_smiles, commercial_dataset))
    
    return output_hits


def target_bond_breaker(target_molecule_smiles, target_substructure_smarts, substructure_atom_map, \
                        disconnection_map, commercial_dataset, building_blocks_a, building_blocks_b, \
                        building_blocks_a_labels, building_blocks_b_labels, output_path, debug):
    # 1. find all target substructures in target
    # 2. make copy of target, number atoms accordingly
    # 3. apply reverse transformation to target
    # 4. save breakdown in text file
    # 5. repeat 2-4 for each target bond found in target
    # return None if no target substructures found in target
    # INPUT: target structure, numbered target substructure, bond to be made
    
    start_time = time.time()
    
    output_file = {
        "input_target_molecule_smiles": [],
        "mapped_input_target_molecule_smiles": [],
        "target_substructure_smarts": [],
        "substructure_atom_map": [],
        "disconnection_map": [],
        "synthon_a": [],
        "synthon_a_label": [],
        "synthon_b": [],
        "synthon_b_label": [],
        "synthon_a_building_block_class": [],
        "synthon_a_building_block_label": [],
        "synthon_b_building_block_class": [],
        "synthon_b_building_block_label": [],
        "synthon_a_building_block_smiles": [],
        "synthon_a_building_block_cas": [],
        "synthon_a_building_block_price": [],
        "synthon_a_building_block_material": [],
        "synthon_b_building_block_smiles": [],
        "synthon_b_building_block_cas": [],
        "synthon_b_building_block_price": [],
        "synthon_b_building_block_material": []
    }
    
    target_mol, message = check_input_target_molecule(target_molecule_smiles)
    if not target_mol: 
        if debug:
            print(message)
        df = pd.DataFrame(output_file)
        df.to_excel(output_path)
        print(f"Execution failed.")
        return None
    
    substructure_mol, message = check_input_target_substructure(target_substructure_smarts)
    if not substructure_mol: 
        if debug:
            print(message)
        df = pd.DataFrame(output_file)
        df.to_excel(output_path)
        print(f"Execution failed.")
        return None

    hits, message = check_number_of_substructure_hits(target_mol, substructure_mol)
    if not hits:
        if debug:
            print(message)
        df = pd.DataFrame(output_file)
        df.to_excel(output_path)
        print(f"Execution failed.")
        return None

    mapped_mols = []
    mapped_mol_smiles_list = []
    for hit in hits:
        mapped_mol, mapped_mol_smiles = copy_and_map_atoms_to_target(target_mol, hit, substructure_atom_map)
        mapped_mols.append(mapped_mol)
        mapped_mol_smiles_list.append(mapped_mol_smiles)
        

    for mapped_sm in mapped_mol_smiles_list:
        # try:
        syn_a, syn_b = disconnect_synthon(mapped_sm, disconnection_map)
        # except:
            # continue
        syn_a_availability = check_synthon_for_commercial_analogs(syn_a, 1, commercial_dataset, building_blocks_a, 3)
        syn_b_availability = check_synthon_for_commercial_analogs(syn_b, 2, commercial_dataset, building_blocks_b, 4)

        for i,sa in enumerate(syn_a_availability):
            for i2,sb in enumerate(syn_b_availability):
                if sa is not None and sb is not None:
                    output_file["input_target_molecule_smiles"].append(target_molecule_smiles)
                    output_file["mapped_input_target_molecule_smiles"].append(mapped_sm)
                    output_file["target_substructure_smarts"].append(target_substructure_smarts)
                    output_file["substructure_atom_map"].append(substructure_atom_map)
                    output_file["disconnection_map"].append(disconnection_map)
                    output_file["synthon_a"].append(syn_a)
                    output_file["synthon_b"].append(syn_b)
                    output_file["synthon_a_label"].append(1)
                    output_file["synthon_b_label"].append(2)
                    output_file["synthon_a_building_block_smiles"].append(sa["sanitized_smiles"])
                    output_file["synthon_a_building_block_cas"].append(sa["cas"])
                    output_file["synthon_a_building_block_price"].append(sa["price"])
                    output_file["synthon_a_building_block_material"].append(sa["material"])
                    output_file["synthon_b_building_block_smiles"].append(sb["sanitized_smiles"])
                    output_file["synthon_b_building_block_cas"].append(sb["cas"])
                    output_file["synthon_b_building_block_price"].append(sb["price"])
                    output_file["synthon_b_building_block_material"].append(sb["material"])
                    output_file["synthon_a_building_block_class"].append(building_blocks_a_labels[i])
                    output_file["synthon_b_building_block_class"].append(building_blocks_b_labels[i2])
                    output_file["synthon_a_building_block_label"].append(3)
                    output_file["synthon_b_building_block_label"].append(4)
    df = pd.DataFrame(output_file)
    df.to_excel(output_path)
    
    end_time = time.time()
    print(f"Found {len(df)} hits for input {target_molecule_smiles} in {str(end_time - start_time)} seconds")


def run(data):
    commercial_data = pd.read_excel(data["commercial_data_path"])

    target_bond_breaker(data["target_molecule"],
                        data["target_substructure"], 
                        data["substructure_atom_map"], 
                        data["disconnection_map"],
                        commercial_data, 
                        data["building_blocks_a"], 
                        data["building_blocks_b"],
                        data["building_blocks_labels_a"], 
                        data["building_blocks_labels_b"], 
                        data["output_path"], 
                        data["debug_flag"])

def main():
    if len(sys.argv) != 2:
        print("Usage: {} <json_file>".format(sys.argv[0]))
        sys.exit(1)

    input_file = sys.argv[1]
    with open(input_file) as f:
        data = json.load(f)
    
    commercial_data = pd.read_excel(data["commercial_data_path"])

    target_bond_breaker(data["target_molecule"],
                        data["target_substructure"], 
                        data["substructure_atom_map"], 
                        data["disconnection_map"],
                        commercial_data, 
                        data["building_blocks_a"], 
                        data["building_blocks_b"],
                        data["building_blocks_labels_a"], 
                        data["building_blocks_labels_b"], 
                        data["output_path"], 
                        data["debug_flag"])
    
if __name__ == "__main__":
    main()
    