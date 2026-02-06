"""Fill cavities with water molecules using Monte Carlo sampling and MMFF94 optimization."""

import os
from typing import List, Dict, Any, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Van der Waals radii for clash detection (in Angstroms)
VDW_RADII = {
    'H': 1.2,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'S': 1.8,
    'P': 1.8,
}

# Minimum distance between water oxygens (in Angstroms)
MIN_WATER_WATER_DISTANCE = 2.8

# Minimum distance between water and protein atoms (in Angstroms)
MIN_WATER_PROTEIN_DISTANCE = 2.4


def fill_cavities_with_water(
    protein_file: str,
    selected_cavities: List[Dict[str, Any]],
    cavity_data: Any,
    output_dir: str,
    waters_per_cavity: Dict[int, int] = None,
) -> str:
    """
    Fill selected cavities with water molecules using Monte Carlo sampling and MMFF94 optimization.
    
    Args:
        protein_file: Path to the protein PDB file
        selected_cavities: List of selected cavity dictionaries
        cavity_data: Cavity data from pyKVFinder
        output_dir: Directory to save output files
        waters_per_cavity: Dictionary mapping cavity ID to number of waters to place
        
    Returns:
        Path to the output PDB file with filled cavities
    """
    from cavefiller.cavity_finder import get_cavity_grid_points
    
    # Output file for the filled protein
    output_file = os.path.join(output_dir, "protein_filled.pdb")
    
    # Read protein structure
    protein_atoms = read_protein_atoms(protein_file)
    
    # Fill each cavity with water
    all_water_molecules = []
    
    for cavity in selected_cavities:
        cavity_id = cavity["id"]
        volume = cavity["volume"]
        
        # Determine number of waters to place
        if waters_per_cavity and cavity_id in waters_per_cavity:
            n_waters = waters_per_cavity[cavity_id]
        else:
            # Default: estimate based on volume
            n_waters = max(1, int(volume / 30))  # ~30 Å³ per water
        
        print(f"Placing {n_waters} water molecules in cavity {cavity_id} (volume: {volume:.2f} Ų)...")
        
        # Get cavity grid points
        cavity_points = get_cavity_grid_points(cavity_data, cavity_id)
        
        if len(cavity_points) == 0:
            print(f"  Warning: No grid points found for cavity {cavity_id}")
            continue
        
        # Place waters using Monte Carlo sampling
        water_positions = monte_carlo_water_placement(
            cavity_points,
            protein_atoms,
            n_waters,
            max_attempts=1000
        )
        
        if len(water_positions) > 0:
            # Optimize water positions using MMFF94
            optimized_waters = optimize_waters_mmff94(
                water_positions,
                protein_atoms,
                cavity_points
            )
            
            all_water_molecules.extend(optimized_waters)
            print(f"  ✅ Successfully placed and optimized {len(optimized_waters)} waters")
        else:
            print(f"  ⚠️  Warning: Could not place any waters in cavity {cavity_id}")
    
    # Write output PDB file
    write_output_pdb(protein_file, all_water_molecules, output_file)
    
    print(f"\nTotal waters added: {len(all_water_molecules)}")
    
    return output_file


def read_protein_atoms(protein_file: str) -> List[Tuple[str, np.ndarray]]:
    """
    Read protein atoms from PDB file.
    
    Args:
        protein_file: Path to PDB file
        
    Returns:
        List of tuples (element, coordinates)
    """
    atoms = []
    
    with open(protein_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                # Parse PDB line
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip()
                    
                    # If element not specified, guess from atom name
                    if not element:
                        atom_name = line[12:16].strip()
                        element = atom_name[0]
                    
                    atoms.append((element, np.array([x, y, z])))
                except (ValueError, IndexError):
                    continue
    
    return atoms


def check_clash(position: np.ndarray, existing_atoms: List[Tuple[str, np.ndarray]], 
                existing_waters: List[np.ndarray]) -> bool:
    """
    Check if a water position clashes with protein atoms or other waters.
    
    Args:
        position: Proposed water oxygen position (x, y, z)
        existing_atoms: List of (element, coordinates) for protein atoms
        existing_waters: List of existing water oxygen positions
        
    Returns:
        True if there is a clash, False otherwise
    """
    # Check clash with protein atoms
    for element, atom_coords in existing_atoms:
        distance = np.linalg.norm(position - atom_coords)
        
        # Get VDW radius for the atom (default to 1.7 if unknown)
        vdw_radius = VDW_RADII.get(element, 1.7)
        
        # Minimum allowed distance is sum of VDW radii with some tolerance
        min_distance = vdw_radius + VDW_RADII['O'] - 0.5  # 0.5 Å tolerance
        
        if distance < max(min_distance, MIN_WATER_PROTEIN_DISTANCE):
            return True
    
    # Check clash with other waters
    for water_pos in existing_waters:
        distance = np.linalg.norm(position - water_pos)
        if distance < MIN_WATER_WATER_DISTANCE:
            return True
    
    return False


def monte_carlo_water_placement(
    cavity_points: np.ndarray,
    protein_atoms: List[Tuple[str, np.ndarray]],
    n_waters: int,
    max_attempts: int = 1000
) -> List[np.ndarray]:
    """
    Place water molecules in cavity using Monte Carlo sampling with clash detection.
    
    Args:
        cavity_points: Array of (x, y, z) coordinates defining the cavity
        protein_atoms: List of (element, coordinates) for protein atoms
        n_waters: Number of water molecules to place
        max_attempts: Maximum attempts per water molecule
        
    Returns:
        List of water oxygen positions (x, y, z)
    """
    placed_waters = []
    
    # Get cavity bounding box
    min_coords = cavity_points.min(axis=0)
    max_coords = cavity_points.max(axis=0)
    
    for i in range(n_waters):
        placed = False
        
        for attempt in range(max_attempts):
            # Random position within cavity bounding box
            position = np.random.uniform(min_coords, max_coords)
            
            # Check if position is within cavity (close to a cavity grid point)
            distances_to_cavity = np.linalg.norm(cavity_points - position, axis=1)
            min_dist_to_cavity = distances_to_cavity.min()
            
            # Position should be within 1.5 Å of a cavity grid point
            if min_dist_to_cavity > 1.5:
                continue
            
            # Check for clashes
            if not check_clash(position, protein_atoms, placed_waters):
                placed_waters.append(position)
                placed = True
                break
        
        if not placed:
            # Could not place this water after max_attempts
            print(f"    Warning: Could only place {len(placed_waters)}/{n_waters} waters")
            break
    
    return placed_waters


def optimize_waters_mmff94(
    water_positions: List[np.ndarray],
    protein_atoms: List[Tuple[str, np.ndarray]],
    cavity_points: np.ndarray,
    max_iterations: int = 200
) -> List[np.ndarray]:
    """
    Optimize water positions using MMFF94 force field.
    
    Args:
        water_positions: Initial water oxygen positions
        protein_atoms: Protein atom coordinates (for reference/constraints)
        cavity_points: Cavity grid points (to ensure waters stay in cavity)
        max_iterations: Maximum optimization iterations
        
    Returns:
        Optimized water oxygen positions
    """
    if len(water_positions) == 0:
        return []
    
    try:
        # Create an rdkit molecule with all water molecules
        mol = Chem.RWMol()
        
        # Add water molecules (H-O-H)
        water_oxygen_indices = []
        for pos in water_positions:
            # Add oxygen
            o_idx = mol.AddAtom(Chem.Atom(8))  # Oxygen
            water_oxygen_indices.append(o_idx)
            
            # Add hydrogens in tetrahedral geometry
            h1_idx = mol.AddAtom(Chem.Atom(1))  # Hydrogen
            h2_idx = mol.AddAtom(Chem.Atom(1))  # Hydrogen
            
            # Add bonds
            mol.AddBond(o_idx, h1_idx, Chem.BondType.SINGLE)
            mol.AddBond(o_idx, h2_idx, Chem.BondType.SINGLE)
        
        # Convert to molecule
        mol = mol.GetMol()
        
        # Sanitize molecule to calculate implicit valences
        Chem.SanitizeMol(mol)
        
        # Create a conformer with the initial positions
        conf = Chem.Conformer(mol.GetNumAtoms())
        
        atom_idx = 0
        for pos in water_positions:
            # Oxygen position
            conf.SetAtomPosition(atom_idx, (float(pos[0]), float(pos[1]), float(pos[2])))
            atom_idx += 1
            
            # Hydrogen positions (tetrahedral geometry around oxygen)
            # H1 at ~0.96 Å from O
            h1_pos = pos + np.array([0.757, 0.586, 0.0])
            conf.SetAtomPosition(atom_idx, (float(h1_pos[0]), float(h1_pos[1]), float(h1_pos[2])))
            atom_idx += 1
            
            # H2 at ~0.96 Å from O
            h2_pos = pos + np.array([-0.757, 0.586, 0.0])
            conf.SetAtomPosition(atom_idx, (float(h2_pos[0]), float(h2_pos[1]), float(h2_pos[2])))
            atom_idx += 1
        
        mol.AddConformer(conf)
        
        # Set up MMFF94 force field
        # Generate MMFF94 properties
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
        
        # Create force field
        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=0)
        
        # Optimize
        ff.Initialize()
        converged = ff.Minimize(maxIts=max_iterations)
        
        # Extract optimized oxygen positions
        optimized_positions = []
        conf = mol.GetConformer(0)
        
        for o_idx in water_oxygen_indices:
            pos = conf.GetAtomPosition(o_idx)
            optimized_positions.append(np.array([pos.x, pos.y, pos.z]))
        
        # Filter waters that moved too far from cavity
        filtered_positions = []
        for pos in optimized_positions:
            distances_to_cavity = np.linalg.norm(cavity_points - pos, axis=1)
            min_dist_to_cavity = distances_to_cavity.min()
            
            # Keep water if it's still reasonably close to cavity
            if min_dist_to_cavity < 2.0:
                filtered_positions.append(pos)
        
        return filtered_positions
        
    except Exception as e:
        print(f"    Warning: MMFF94 optimization failed: {e}")
        print(f"    Returning unoptimized positions")
        return water_positions


def write_output_pdb(
    protein_file: str,
    water_positions: List[np.ndarray],
    output_file: str
) -> None:
    """
    Write protein and water molecules to output PDB file.
    
    Args:
        protein_file: Path to input protein PDB file
        water_positions: List of water oxygen positions
        output_file: Path to output PDB file
    """
    # Read protein structure
    with open(protein_file, 'r') as f:
        protein_lines = f.readlines()
    
    # Write protein atoms (exclude END line)
    with open(output_file, 'w') as f:
        for line in protein_lines:
            if not line.startswith('END'):
                f.write(line)
        
        # Add water molecules
        for i, pos in enumerate(water_positions, start=1):
            # Write oxygen
            water_line = (
                f"HETATM{i:5d}  O   HOH W{i:4d}    "
                f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00           O\n"
            )
            f.write(water_line)
        
        f.write('END\n')
