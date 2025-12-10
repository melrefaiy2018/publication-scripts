"""
PARALLEL-SAFE VERSION: Includes seed in Hamiltonian filenames to avoid race conditions
when running multiple jobs in parallel.

Each job creates its own unique Hamiltonian file based on seed.
"""

import warnings
import os
import copy
import numpy as np
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio.PDB import PDBParser, Chain
from pymembrane.structure.atomic_pigment import ChlorophyllAtomic
from pymembrane.structure.atomic_protein import ElectrostaticProteinAtomic
from pymembrane.parameters.thylakoid.plants.hamiltonian.pigment_renger_hfcis import (
    coupling_data_by_type_mcce as coupling_renger_by_type_mcce,
    q0011_charges_by_type_mcce as q0011_renger_by_type_mcce
)
from pymembrane.parameters.thylakoid.cyano.hamiltonian.pigment_renger_lineshape_isia import kubo_renger_cla_reduce
from pymembrane.parameters.thylakoid.cyano.hamiltonian.pigment_novo_lineshape import kubo_novo_cla_ct


def validate_pdb_structure(pdb_file, verbose=True):
    """
    Validate PDB file has all required pigments before processing.
    """
    if not os.path.exists(pdb_file):
        raise ValueError(f"PDB file not found: {pdb_file}")
    
    pdb_parser = PDBParser(QUIET=True)
    try:
        pdb_atomic = pdb_parser.get_structure('Isia', pdb_file)
        model = pdb_atomic[0]
    except Exception as e:
        raise ValueError(f"Failed to parse PDB file {pdb_file}: {e}")
    
    # Check Chain A exists
    if 'A' not in model.child_dict:
        available = list(model.child_dict.keys())
        raise ValueError(f"Chain A not found in {os.path.basename(pdb_file)}. Available chains: {available}")
    
    chain_A = model['A']
    
    # Extract CLA residues
    cla_residues = {}
    for residue in chain_A:
        res_id = residue.get_id()
        hetero_flag = res_id[0].strip()
        res_num = res_id[1]
        
        if hetero_flag in ['CLA', 'CL', 'H_CLA', ''] and residue.resname in ['CLA', 'CL']:
            cla_residues[res_num] = hetero_flag if hetero_flag else residue.resname
    
    if not cla_residues:
        raise ValueError(f"No CLA residues found in Chain A of {os.path.basename(pdb_file)}")
    
    cla_numbers = sorted(cla_residues.keys())
    
    if verbose:
        print(f"✓ PDB Validation: Chain A found with CLA residues {cla_numbers}")
    
    return cla_residues


def construct_monomer(pdb_file, dir_path=None, e0_a=14900, e0_b=0, N_ens=100, temp=300, 
                     dielectric_tresp=1, dielectric_dipole=1, dielectric_cdc=2, seed=None):
    """
    Construct an IsiA monomer without charge transfer state.
    
    PARALLEL-SAFE: Includes seed in filename to avoid conflicts
    
    Args:
        pdb_file (str): Path to the PDB file
        dir_path (str, optional): Directory path for output files
        e0_a (float): Site energy parameter a
        e0_b (float): Site energy parameter b
        N_ens (int): Number of ensembles
        temp (float): Temperature in Kelvin
        dielectric_tresp (float): Dielectric constant for TRESP
        dielectric_dipole (float): Dielectric constant for dipole
        dielectric_cdc (float): Dielectric constant for CDC
        seed (int, optional): Random seed - if provided, included in Hamiltonian filename
                             for parallel processing safety
    
    Returns:
        ElectrostaticProteinAtomic: Configured IsiA atomic object
    """
    
    # Set directory path
    if dir_path is None:
        dir_path = os.getcwd()
    
    # Validate PDB structure first
    try:
        validate_pdb_structure(pdb_file, verbose=True)
    except ValueError as e:
        raise ValueError(f"PDB validation failed: {e}")
    
    # Load PDB structure
    try:
        pdb_parser = PDBParser(QUIET=True)
        pdb_atomic = pdb_parser.get_structure('Isia', pdb_file)
        model = pdb_atomic[0]
    except Exception as e:
        raise ValueError(f"Failed to load PDB file '{pdb_file}': {e}")

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
        
        # Replace NaN with zero in dict
        for inner_dict in df_hamiltonian.values():
            for key, value in inner_dict.items():
                if isinstance(value, float) and np.isnan(value):
                    inner_dict[key] = 0
        
        # Save Hamiltonian - PARALLEL-SAFE: include seed in filename
        os.makedirs(os.path.join(dir_path, 'Hamiltonian'), exist_ok=True)
        
        # Add seed to filename if provided (for parallel safety)
        seed_str = f"_seed{seed}" if seed is not None else ""
        path_hamiltonian = os.path.join(
            dir_path, 
            f'Hamiltonian/fit_IsiA_chain_A_dielc_{IsiA_atomic.dielectric_cdc}_a{e0_a}_b{e0_b}{seed_str}.csv'
        )
        
        IsiA_atomic._save_hamiltonian_to_csv(df_hamiltonian, path_hamiltonian)
        
        # Verify file exists and has content
        if not os.path.exists(path_hamiltonian):
            raise ValueError(f"Hamiltonian file was not created: {path_hamiltonian}")
        file_size = os.path.getsize(path_hamiltonian)
        if file_size < 500:
            raise ValueError(f"Hamiltonian file suspiciously small ({file_size} bytes)")
        
        # Load Hamiltonian with error handling
        try:
            IsiA_atomic.load_hamiltonian(path_hamiltonian)
        except IndexError as e:
            raise ValueError(
                f"Hamiltonian loading failed - pigment name mismatch. "
                f"This means pigment names in the CSV file don't match the pigments "
                f"currently detected in the protein structure. Original error: {e}"
            )
        
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
        IsiA_atomic.define_domains(list_pigments_by_domain)

        return IsiA_atomic
    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Failed to construct IsiA atomic object: {e}. "
                        f"This may be due to: (1) missing Chain A in PDB file, "
                        f"(2) improper pigment preparation/parameters, (3) Hamiltonian construction issues, "
                        f"or (4) invalid dielectric constants.")
    

def construct_monomer_with_CT(pdb_file, dir_path=None, e0_a=14900, ct_energy=14580, e0_b=0, N_ens=100, temp=300,
                           dielectric_tresp=1, dielectric_dipole=1.5, dielectric_cdc=2, coupling_ct=50,
                           ct_pigment='502', seed=None):
    """
    Construct an IsiA monomer with charge transfer state.
    
    PARALLEL-SAFE: Includes seed in filename to avoid conflicts
    
    Args:
        pdb_file (str): Path to the PDB file
        dir_path (str, optional): Directory path for output files
        e0_a (float): Site energy parameter a
        ct_energy (float): CT state energy
        e0_b (float): Site energy parameter b
        N_ens (int): Number of ensembles
        temp (float): Temperature in Kelvin
        dielectric_tresp (float): Dielectric constant for TRESP
        dielectric_dipole (float): Dielectric constant for dipole
        dielectric_cdc (float): Dielectric constant for CDC
        coupling_ct (float): CT coupling value
        ct_pigment (str): CT pigment identifier
        seed (int, optional): Random seed - if provided, included in Hamiltonian filename
                             for parallel processing safety
    
    Returns:
        ElectrostaticProteinAtomic: Configured IsiA atomic object
    """
    
    # Set directory path
    if dir_path is None:
        dir_path = os.getcwd()
    
    # Validate PDB structure first
    try:
        validate_pdb_structure(pdb_file, verbose=True)
    except ValueError as e:
        raise ValueError(f"PDB validation failed: {e}")
    
    # Load PDB structure
    try:
        pdb_parser = PDBParser(QUIET=True)
        pdb_atomic = pdb_parser.get_structure('Isia', pdb_file)
        model = pdb_atomic[0]
    except Exception as e:
        raise ValueError(f"Failed to load PDB file '{pdb_file}': {e}")

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
            raise ValueError(f"Residue {residue_id} not found in Chain A. The PDB file must contain "
                           f"a CLA residue (Chlorophyll A) with residue number 502 in Chain A.")
    except Exception as e:
        raise ValueError(f"Failed to process Chain A and create charge transfer state: {e}")

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
        
        # Save Hamiltonian - PARALLEL-SAFE: include seed in filename
        os.makedirs(os.path.join(dir_path, 'Hamiltonian'), exist_ok=True)
        
        # Add seed to filename if provided (for parallel safety)
        seed_str = f"_seed{seed}" if seed is not None else ""
        path_hamiltonian = os.path.join(
            dir_path, 
            f'Hamiltonian/fit_CT_IsiA_chain_A_dielc_{IsiA_atomic_CT.dielectric_cdc}_a{e0_a}_b{e0_b}_coupling_{coupling_ct}{seed_str}.csv'
        )
        
        IsiA_atomic_CT._save_hamiltonian_to_csv(df_hamiltonian, path_hamiltonian)
        
        # Verify file exists and has content
        if not os.path.exists(path_hamiltonian):
            raise ValueError(f"Hamiltonian file was not created: {path_hamiltonian}")
        file_size = os.path.getsize(path_hamiltonian)
        if file_size < 500:
            raise ValueError(f"Hamiltonian file suspiciously small ({file_size} bytes)")
        
        # Load Hamiltonian with error handling
        try:
            IsiA_atomic_CT.load_hamiltonian(path_hamiltonian)
        except IndexError as e:
            raise ValueError(
                f"CT Hamiltonian loading failed - pigment name mismatch. "
                f"This means pigment names in the CSV file don't match the pigments "
                f"currently detected in the protein structure (Chains A + Z). Original error: {e}"
            )
        
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
        
    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"Failed to construct IsiA atomic object with charge transfer state: {e}. "
                        f"Possible causes: (1) Chain Z CT pigment preparation failed, "
                        f"(2) Hamiltonian construction or coupling parameters invalid, "
                        f"(3) incompatible lineshape function (kubo_novo_cla_ct), "
                        f"(4) ct_pigment identifier '{ct_pigment}' does not exist, "
                        f"or (5) invalid dielectric constants or other parameters.")
