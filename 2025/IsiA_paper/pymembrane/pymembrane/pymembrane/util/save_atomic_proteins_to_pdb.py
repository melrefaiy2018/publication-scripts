"""
Save two proteins to a single PDB file for visualization in Chimera/PyMOL
This allows you to see the true 3D structure and Z-axis separation
Uses ExtendedPDBWriter to preserve charge information
"""
from Bio.PDB import PDBIO, Structure, Model, Chain
from Bio.PDB.Atom import Atom
import numpy as np
import sys
from pathlib import Path

# Add pymembrane to path
pymembrane_path = Path(__file__).parent.parent.parent.parent.parent.parent.parent / "pymembrane"
if str(pymembrane_path) not in sys.path:
    sys.path.insert(0, str(pymembrane_path))

from pymembrane.util.writeExtendedPDB import ExtendedPDBWriter

def save_proteins_to_pdb(proteins, output_file):
    """
    Save a list of CGProtein objects to a single PDB file using ExtendedPDBWriter
    This preserves charge information and creates a proper extended PDB
    Each protein is assigned to a different chain (A, B, C, ...)

    Parameters:
    -----------
    proteins : list of CGProtein
        List of protein objects to save
    output_file : str or Path
        Output PDB filename
    """

    # Remove the file if it exists (ExtendedPDBWriter appends)
    output_file = Path(output_file)
    if output_file.exists():
        output_file.unlink()

    # Create writer
    writer = ExtendedPDBWriter(str(output_file))

    atom_counter = 1
    total_atoms = 0
    protein_stats = []
    
    # Track chain IDs used across all proteins to avoid collisions
    used_chain_ids = set()

    # Loop through each protein and preserve its original chain structure
    for prot_idx, protein in enumerate(proteins):
        # Get atoms and coordinates
        prot_atoms = list(protein.atomic.get_atoms())
        prot_coords = protein.R2_atomic_xyz

        # Count chains in this protein
        chain_count = {}
        for atom in prot_atoms:
            orig_chain_id = atom.get_parent().get_parent().id
            chain_count[orig_chain_id] = chain_count.get(orig_chain_id, 0) + 1
        
        print(f"Adding {protein.name} ({len(prot_atoms)} atoms) with {len(chain_count)} chain(s): {list(chain_count.keys())}...")

        # Write protein preserving original chain IDs
        for atom, coord in zip(prot_atoms, prot_coords):
            residue = atom.get_parent()
            orig_chain = residue.get_parent()
            orig_chain_id = orig_chain.id
            
            # Use original chain ID (already unique per protein in most cases)
            # If there's a conflict, we'll handle it
            chain_id = orig_chain_id
            if chain_id in used_chain_ids:
                # Find an alternative chain ID
                import string
                available_ids = (list(string.ascii_uppercase) + 
                               list(string.digits) + 
                               list(string.ascii_lowercase))
                for alt_id in available_ids:
                    if alt_id not in used_chain_ids:
                        chain_id = alt_id
                        break
            
            used_chain_ids.add(chain_id)

            # Handle magnesium element name
            if atom.get_name() == 'MG':
                atom_element = 'Mg'
            else:
                atom_element = atom.element

            # Get charge (default to 0 if not available)
            charge = atom.get_bfactor() if hasattr(atom, 'get_bfactor') else 0.0
            if hasattr(protein.atomic, '_dict_atom_charges') and protein.atomic._dict_atom_charges:
                atom_key = f"{orig_chain.id}_{residue.resname}_{residue.id[1]}_{atom.get_name()}"
                charge = protein.atomic._dict_atom_charges.get(atom_key, 0.0)

            atom_info = {
                "atom_record": "ATOM",
                "atom_index": atom_counter,
                "atom_name": atom.get_name(),
                "residue_name": residue.resname,
                "chain_id": chain_id,
                "residue_index": str(residue.id[1]),
                "x": float(coord[0]),
                "y": float(coord[1]),
                "z": float(coord[2]),
                "element_name": atom_element,
                "charge": charge,
            }
            writer.write_pdb_line(**atom_info)
            atom_counter += 1

        total_atoms += len(prot_atoms)
        protein_stats.append((list(chain_count.keys()), protein.name, len(prot_atoms)))

    # Print summary
    print(f"\nExtended PDB file saved: {output_file}")
    print(f"  Total proteins: {len(proteins)}")
    print(f"  Total atoms: {total_atoms}")
    print(f"\nProtein chain assignments:")
    for chain_ids, prot_name, n_atoms in protein_stats:
        chains_str = ", ".join(chain_ids)
        print(f"  {prot_name}: Chains {chains_str} ({n_atoms} atoms)")


def save_two_proteins_to_pdb(prot1, prot2, output_file):
    """
    Convenience function to save two proteins to a PDB file
    This is a wrapper around save_proteins_to_pdb for backward compatibility

    Parameters:
    -----------
    prot1 : CGProtein
        First protein
    prot2 : CGProtein
        Second protein
    output_file : str or Path
        Output PDB filename
    """
    save_proteins_to_pdb([prot1, prot2], output_file)

    # Print helpful information for visualization
    print(f"\n{'='*70}")
    print("HOW TO VISUALIZE IN CHIMERA/PYMOL")
    print(f"{'='*70}")

    # Calculate centers and distances
    prot1_center = prot1.location
    prot2_center = prot2.location
    distance = np.linalg.norm(prot1_center - prot2_center)

    # Calculate Z-axis statistics
    prot1_z_coords = prot1.R2_atomic_xyz[:, 2]
    prot2_z_coords = prot2.R2_atomic_xyz[:, 2]

    prot1_z_range = (np.min(prot1_z_coords), np.max(prot1_z_coords))
    prot2_z_range = (np.min(prot2_z_coords), np.max(prot2_z_coords))

    z_separation = np.mean(prot1_z_coords) - np.mean(prot2_z_coords)

    print(f"\nProtein Information:")
    print(f"  Protein 1 ({prot1.name}):")
    print(f"    Center: {prot1_center}")
    print(f"    Z range: {prot1_z_range[0]:.2f} to {prot1_z_range[1]:.2f} Å")
    print(f"    Mean Z: {np.mean(prot1_z_coords):.2f} Å")
    
    # Show chains in prot1
    prot1_chains = set(atom.get_parent().get_parent().id for atom in prot1.atomic.get_atoms())
    print(f"    Chains: {', '.join(sorted(prot1_chains))}")

    print(f"\n  Protein 2 ({prot2.name}):")
    print(f"    Center: {prot2_center}")
    print(f"    Z range: {prot2_z_range[0]:.2f} to {prot2_z_range[1]:.2f} Å")
    print(f"    Mean Z: {np.mean(prot2_z_coords):.2f} Å")
    
    # Show chains in prot2
    prot2_chains = set(atom.get_parent().get_parent().id for atom in prot2.atomic.get_atoms())
    print(f"    Chains: {', '.join(sorted(prot2_chains))}")

    print(f"\n  Distance between centers: {distance:.2f} Å")
    print(f"  Z-axis separation (mean): {abs(z_separation):.2f} Å")

    print(f"\n{'='*70}")
    print("CHIMERA COMMANDS:")
    print(f"{'='*70}")
    
    all_chains = sorted(prot1_chains | prot2_chains)
    prot1_chain_str = " or :." + " or :.".join(sorted(prot1_chains))
    prot2_chain_str = " or :." + " or :.".join(sorted(prot2_chains))
    
    print(f"""
# Open the file
open {output_file}

# Color proteins differently (preserving original chain structure)
color red {prot1_chain_str}
color blue {prot2_chain_str}

# Show as spheres
represent sphere

# View from top (looking down Z-axis)
view orient
turn z 0

# Measure distance between protein centers (pick representative atoms)
# You can manually measure using the tape measure tool

# To see the 2D vs 3D effect:
# 1. Top view (looks overlapped): turn z 0
# 2. Side view (shows Z-separation): turn x 90

# All chains present: {', '.join(all_chains)}
# Protein 1 ({prot1.name}): {', '.join(sorted(prot1_chains))}
# Protein 2 ({prot2.name}): {', '.join(sorted(prot2_chains))}
    """)

    print(f"\n{'='*70}")
    print("PYMOL COMMANDS:")
    print(f"{'='*70}")
    
    prot1_chain_sel = "+".join(sorted(prot1_chains))
    prot2_chain_sel = "+".join(sorted(prot2_chains))
    
    print(f"""
# Load the file
load {output_file}

# Color proteins (preserving original chains)
color red, chain {prot1_chain_sel}
color blue, chain {prot2_chain_sel}

# Show as spheres
show spheres

# Top view (XY plane)
orient
turn z, 0

# Side view to see Z-separation
# turn x, 90

# Measure distance between protein centers
# Use measurement tool or distance command

# All chains: {', '.join(all_chains)}
# Protein 1 ({prot1.name}): chains {prot1_chain_sel}
# Protein 2 ({prot2.name}): chains {prot2_chain_sel}
    """)

    print(f"\n{'='*70}")
    print("WHAT TO LOOK FOR:")
    print(f"{'='*70}")
    print(f"""
1. TOP VIEW (looking down Z-axis):
   → Chains will appear close together in XY plane
   → This mimics your 2D matplotlib plots
   → May look like they're overlapping!

2. SIDE VIEW (rotate 90° around X):
   → You'll see the Z-axis separation clearly
   → Chain A (red) at Z ≈ {np.mean(prot1_z_coords):.1f} Å
   → Chain B (blue) at Z ≈ {np.mean(prot2_z_coords):.1f} Å
   → Vertical gap ≈ {abs(z_separation):.1f} Å

3. 3D ROTATION:
   → Rotate the structure freely
   → Confirm there's NO atomic overlap
   → See how 2D projection can be misleading!
    """)



# # ======================================================================
# # SAVE PROTEINS TO PDB FOR CHIMERA/PYMOL VISUALIZATION
# # ======================================================================

# from save_proteins_to_pdb import save_two_proteins_to_pdb, save_protein_with_highlighted_pigments

# print(f"\n{'='*70}")
# print("SAVING PROTEINS TO PDB FILE FOR 3D VISUALIZATION")
# print(f"{'='*70}\n")

# # Save the two proteins to a PDB file
# pdb_output = dir_save / f'proteins_{prot1.name}_{prot2.name}_seed{seed}.pdb'
# save_two_proteins_to_pdb(prot1, prot2, pdb_output)

# # Also save with highlighted pigments
# pdb_highlighted = dir_save / f'proteins_highlighted_{prot1.name}_{prot2.name}_seed{seed}.pdb'
# save_protein_with_highlighted_pigments(
#     prot1, prot2,
#     'c2s2_7_5xnl_c2s2_g_CLA_603',  # Pigment from prot1
#     'lhcii_471_1rwt_D_CHL_609',    # Pigment from prot2
#     pdb_highlighted
# )

# print(f"\n{'='*70}")
# print("PDB FILES READY FOR VISUALIZATION!")
# print(f"{'='*70}")
# print(f"\nStandard PDB: {pdb_output}")
# print(f"Highlighted:  {pdb_highlighted}")
# print(f"\nOpen these files in Chimera or PyMOL to see the true 3D structure")
# print(f"and confirm the Z-axis separation that explains the 2D 'overlap'!")
# print(f"{'='*70}\n")
