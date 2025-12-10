import warnings
from Bio import BiopythonWarning, PDB
warnings.simplefilter('ignore', BiopythonWarning)
from pymembrane.structure.atomic_pigment import ChlorophyllAtomic
from pymembrane.parameters.thylakoid.plants.hamiltonian.pigment_renger_hfcis import (
    coupling_data_by_type_mcce as coupling_renger_by_type_mcce,
    q0011_charges_by_type_mcce as q0011_renger_by_type_mcce
)
from pymembrane.parameters.thylakoid.cyano.hamiltonian.pigment_renger_lineshape_isia import kubo_renger_cla_reduce
from pymembrane.parameters.thylakoid.cyano.hamiltonian.pigment_novo_lineshape import kubo_novo_cla_ct
from pymembrane.structure.atomic_protein import ElectrostaticProteinAtomic
from pymembrane.exciton.linear_spectra import LinearSpectra
from pymembrane.util.helper_functions import *
from pymembrane.util.physical_constants import hbar, c
from pymembrane.util.plot_cdc_contribution import plot_residue_contribution_and_spatial_interactions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import copy
from Bio.PDB import PDBParser, Chain


def construct_monomer(pdb_file, dir_path=None, e0_a=14900,  e0_b=0, N_ens=100, temp=300, dielectric_tresp=1, dielectric_dipole=1, dielectric_cdc=2):
    """
    Construct an IsiA monomer with charge transfer state.
    
    Args:
        pdb_file (str): Path to the PDB file
        dir_path (str, optional): Directory path for output files. Defaults to current directory.
        e0_a (float, optional): Site energy parameter a. Defaults to 14900.
        e0_b (float, optional): Site energy parameter b. Defaults to 0.
        N_ens (int, optional): Number of ensembles. Defaults to 100.
        temp (float, optional): Temperature in Kelvin. Defaults to 300.
        dielectric_tresp (float, optional): Dielectric constant for TRESP. Defaults to 1.
        dielectric_dipole (float, optional): Dielectric constant for dipole. Defaults to 1.5.
        dielectric_cdc (float, optional): Dielectric constant for CDC. Defaults to 2.2.
    
    Returns:
        ElectrostaticProteinAtomic: Configured IsiA atomic object
    """
    
    # Set directory path
    if dir_path is None:
        dir_path = os.getcwd()
    
    # Load PDB structure
    try:
        pdb_parser = PDBParser(QUIET=True)
        pdb_atomic = pdb_parser.get_structure('Isia', pdb_file)
        model = pdb_atomic[0]  # First model
    except Exception as e:
        raise ValueError(f"Error loading PDB file: {e}")

    # Initialize IsiA atomic object
    try:
        IsiA_atomic = ElectrostaticProteinAtomic(pdb_atomic, path_extended_pdb=pdb_file, name='Isia')
        
        # Set protein parameters
        IsiA_atomic.fval_tresp = 0.72
        IsiA_atomic.dielectric_tresp = dielectric_tresp
        IsiA_atomic.dielectric_dipole = dielectric_dipole
        IsiA_atomic.dielectric_cdc = dielectric_cdc
        
        # Prepare pigments for Chain A
        IsiA_atomic.prepare_pigments('CLA', ChlorophyllAtomic,
                                list_chain=['A'],
                                lineshape=kubo_renger_cla_reduce,
                                dict_coupling_data=coupling_renger_by_type_mcce['CLA_IPPC'],
                                dict_q0011_charges=q0011_renger_by_type_mcce['CLA'],
                                disorder=70)
        
        
        # Construct Hamiltonian
        df_hamiltonian = IsiA_atomic.construct_hamiltonian(e0_a=e0_a, e0_b=e0_b)
        # df_hamiltonian['Isia_A_CLA_503']['Isia_A_CLA_503'] += 100
        
        # Replace NaN with zero in dict
        for inner_dict in df_hamiltonian.values():
            for key, value in inner_dict.items():
                if isinstance(value, float) and np.isnan(value):
                    inner_dict[key] = 0
        
        # Save Hamiltonian
        os.makedirs(os.path.join(dir_path, 'Hamiltonian'), exist_ok=True)
        path_hamiltonian = os.path.join(dir_path, f'Hamiltonian/fit_IsiA_chain_A_dielc_{IsiA_atomic.dielectric_cdc}_a{e0_a}_b{e0_b}.csv')
        IsiA_atomic._save_hamiltonian_to_csv(df_hamiltonian, path_hamiltonian)
        IsiA_atomic.load_hamiltonian(path_hamiltonian)
        # list_pigments_by_domain = IsiA_atomic.build_hamiltonian_domains('energy', domain_cutoff=0.1)
        # list_pigments_by_domain = [['Isia_A_CLA_501'], ['Isia_A_CLA_503'], ['Isia_A_CLA_504','Isia_A_CLA_502'], 
        #                            ['Isia_A_CLA_519', 'Isia_A_CLA_511', 'Isia_A_CLA_516', 'Isia_A_CLA_510', 'Isia_A_CLA_512', 'Isia_A_CLA_507', 'Isia_A_CLA_509', 'Isia_A_CLA_505', 'Isia_A_CLA_518', 'Isia_A_CLA_506', 'Isia_A_CLA_508'],
        #                            ['Isia_A_CLA_513'], ['Isia_A_CLA_517']]
        # July29, 2025
        list_pigments_by_domain = [['Isia_A_CLA_501'], ['Isia_A_CLA_504', 'Isia_A_CLA_502'], ['Isia_A_CLA_503'], ['Isia_A_CLA_507',
                                    'Isia_A_CLA_510',
                                    'Isia_A_CLA_509',
                                    'Isia_A_CLA_508',
                                    'Isia_A_CLA_505',
                                    'Isia_A_CLA_511'],
                                    ['Isia_A_CLA_506'],
                                    ['Isia_A_CLA_513', 'Isia_A_CLA_512', 'Isia_A_CLA_516'],
                                    ['Isia_A_CLA_517'],
                                    ['Isia_A_CLA_518'],
                                    ['Isia_A_CLA_519']]
        # list_pigments_by_domain = IsiA_atomic.build_hamiltonian_domains('coupling')
        IsiA_atomic.define_domains(list_pigments_by_domain)

        return IsiA_atomic
    except Exception as e:
        raise ValueError(f"Error constructing IsiA atomic object: {e}")
    

def construct_monomer_with_CT(pdb_file, dir_path=None, e0_a=14900, ct_energy=14580, e0_b=0, N_ens=100, temp=300,
                           dielectric_tresp=1, dielectric_dipole=1.5, dielectric_cdc=2, coupling_ct=50,
                           ct_pigment='502'):
    """
    Construct an IsiA monomer with charge transfer state.
    
    Args:
        pdb_file (str): Path to the PDB file
        dir_path (str, optional): Directory path for output files. Defaults to current directory.
        e0_a (float, optional): Site energy parameter a. Defaults to 14900.
        e0_b (float, optional): Site energy parameter b. Defaults to 0.
        N_ens (int, optional): Number of ensembles. Defaults to 100.
        temp (float, optional): Temperature in Kelvin. Defaults to 300.
        dielectric_tresp (float, optional): Dielectric constant for TRESP. Defaults to 1.
        dielectric_dipole (float, optional): Dielectric constant for dipole. Defaults to 1.5.
        dielectric_cdc (float, optional): Dielectric constant for CDC. Defaults to 2.2.
        ct_pigment (str, optional): CT pigment identifier. Defaults to '502'.
    
    Returns:
        ElectrostaticProteinAtomic: Configured IsiA atomic object
    """
    
    # Set directory path
    if dir_path is None:
        dir_path = os.getcwd()
    
    # Load PDB structure
    try:
        pdb_parser = PDBParser(QUIET=True)
        pdb_atomic = pdb_parser.get_structure('Isia', pdb_file)
        model = pdb_atomic[0]  # First model
    except Exception as e:
        raise ValueError(f"Error loading PDB file: {e}")

    # Select Chain A
    try:
        chain_A = model["A"]
        residue_id = ('H_CLA', 502, ' ')
        
        if residue_id in chain_A:
            original_residue = chain_A[residue_id]
            
            # Create a deep copy of the residue
            dummy_residue = copy.deepcopy(original_residue)
            
            # Modify the copied residue properties
            dummy_residue.resname = "CLA"
            dummy_residue.id = ('H_CLA', 600, ' ')
            
            # Ensure Chain Z exists in the model
            if "Z" not in model.child_dict:
                new_chain = Chain.Chain("Z")
                model.add(new_chain)
            else:
                new_chain = model["Z"]
            
            # Translate the copied pigment
            translation_vector = np.array([10.0, 0.0, 0.0])
            for atom in dummy_residue:
                atom.coord += translation_vector
            
            # Add the copied residue to Chain Z
            new_chain.add(dummy_residue)
        else:
            raise ValueError(f"Residue {residue_id} not found in Chain A")
    except Exception as e:
        raise ValueError(f"Error processing Chain A: {e}")

    # Initialize IsiA atomic object
    try:
        IsiA_atomic_CT = ElectrostaticProteinAtomic(pdb_atomic, path_extended_pdb=pdb_file, name='Isia')
        
        # Set protein parameters
        IsiA_atomic_CT.fval_tresp = 0.72
        IsiA_atomic_CT.dielectric_tresp = dielectric_tresp
        IsiA_atomic_CT.dielectric_dipole = dielectric_dipole
        IsiA_atomic_CT.dielectric_cdc = dielectric_cdc
        
        # Prepare pigments for Chain A
        IsiA_atomic_CT.prepare_pigments('CLA', ChlorophyllAtomic,
                                list_chain=['A'],
                                lineshape=kubo_renger_cla_reduce,
                                dict_coupling_data=coupling_renger_by_type_mcce['CLA_IPPC'],
                                dict_q0011_charges=q0011_renger_by_type_mcce['CLA'],
                                disorder=70)
        
        # Prepare pigments for Chain Z
        IsiA_atomic_CT.prepare_pigments('CLA', ChlorophyllAtomic,
                                list_chain=['Z'],
                                lineshape=kubo_novo_cla_ct,
                                dict_coupling_data=coupling_renger_by_type_mcce['CLA_CT'],
                                dict_q0011_charges=q0011_renger_by_type_mcce['CLA_CT'],
                                disorder=70)
        
        # Construct Hamiltonian
        df_hamiltonian = IsiA_atomic_CT.construct_hamiltonian(e0_a=e0_a, e0_b=e0_b)
        df_hamiltonian['Isia_Z_CLA_600']['Isia_Z_CLA_600'] = ct_energy
        df_hamiltonian['Isia_Z_CLA_600'][f'Isia_A_CLA_{ct_pigment}'] = coupling_ct
        df_hamiltonian[f'Isia_A_CLA_{ct_pigment}']['Isia_Z_CLA_600'] = coupling_ct
        
        # Replace NaN with zero in dict
        for inner_dict in df_hamiltonian.values():
            for key, value in inner_dict.items():
                if isinstance(value, float) and np.isnan(value):
                    inner_dict[key] = 0
        
        # Save Hamiltonian
        os.makedirs(os.path.join(dir_path, 'Hamiltonian'), exist_ok=True)
        path_hamiltonian = os.path.join(dir_path, f'Hamiltonian/fit_CT_IsiA_chain_A_dielc_{IsiA_atomic_CT.dielectric_cdc}_a{e0_a}_b{e0_b}_coupling_{coupling_ct}.csv')
        IsiA_atomic_CT._save_hamiltonian_to_csv(df_hamiltonian, path_hamiltonian)
        IsiA_atomic_CT.load_hamiltonian(path_hamiltonian)
        # list_pigments_by_domain = IsiA_atomic.build_hamiltonian_domains('energy', domain_cutoff=0.1)
        list_pigments_by_domain = [['Isia_A_CLA_501'], ['Isia_A_CLA_504', 'Isia_A_CLA_502', 'Isia_Z_CLA_600'], ['Isia_A_CLA_503'], ['Isia_A_CLA_507',
                            'Isia_A_CLA_510',
                            'Isia_A_CLA_509',
                            'Isia_A_CLA_508',
                            'Isia_A_CLA_505',
                            'Isia_A_CLA_511'],
                            ['Isia_A_CLA_506'],
                            ['Isia_A_CLA_513', 'Isia_A_CLA_512', 'Isia_A_CLA_516'],
                            ['Isia_A_CLA_517'],
                            ['Isia_A_CLA_518'],
                            ['Isia_A_CLA_519']]
        IsiA_atomic_CT.define_domains(list_pigments_by_domain)

        return IsiA_atomic_CT
        
    except Exception as e:
        raise ValueError(f"Error constructing IsiA atomic object: {e}")