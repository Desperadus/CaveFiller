"""Fill cavities with water molecules using Packmol."""

import os
import subprocess
import tempfile
from typing import List, Dict, Any
from pathlib import Path
import numpy as np

# Estimated number of water molecules per 1000 cubic angstroms of cavity volume
# This is a rough approximation based on water density
WATERS_PER_1000_CUBIC_ANGSTROMS = 30

# Maximum number of water molecules to place per cavity (for simple placement)
MAX_WATERS_PER_CAVITY = 50


def fill_cavities_with_water(
    protein_file: str,
    selected_cavities: List[Dict[str, Any]],
    cavity_data: Any,
    output_dir: str,
) -> str:
    """
    Fill selected cavities with water molecules using Packmol.
    
    Args:
        protein_file: Path to the protein PDB file
        selected_cavities: List of selected cavity dictionaries
        cavity_data: Cavity data from pyKVFinder
        output_dir: Directory to save output files
        
    Returns:
        Path to the output PDB file with filled cavities
    """
    from cavefiller.cavity_finder import get_cavity_grid_points
    
    # Output file for the filled protein
    output_file = os.path.join(output_dir, "protein_filled.pdb")
    
    # Check if packmol is available
    if not is_packmol_available():
        print("⚠️  Warning: Packmol is not available in PATH.")
        print("Falling back to simple water placement method.")
        return fill_with_simple_placement(
            protein_file, selected_cavities, cavity_data, output_file
        )
    
    # Get bounding boxes for selected cavities
    cavity_regions = []
    for cavity in selected_cavities:
        points = get_cavity_grid_points(cavity_data, cavity["id"])
        if len(points) > 0:
            min_coords = points.min(axis=0)
            max_coords = points.max(axis=0)
            cavity_regions.append({
                "id": cavity["id"],
                "min": min_coords,
                "max": max_coords,
                "center": (min_coords + max_coords) / 2,
                "volume": cavity["volume"],
            })
    
    if not cavity_regions:
        print("Warning: No valid cavity regions found.")
        # Just copy the protein file
        with open(protein_file, 'r') as f_in, open(output_file, 'w') as f_out:
            f_out.write(f_in.read())
        return output_file
    
    # Create water molecule template (single water)
    water_pdb = create_water_template(output_dir)
    
    # Create Packmol input file
    packmol_input = create_packmol_input(
        protein_file, water_pdb, cavity_regions, output_file
    )
    
    # Run Packmol
    run_packmol(packmol_input, output_dir)
    
    return output_file


def is_packmol_available() -> bool:
    """Check if Packmol is available in the system PATH."""
    try:
        subprocess.run(
            ["packmol", "-h"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5,
        )
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        return False


def create_water_template(output_dir: str) -> str:
    """Create a PDB file with a single water molecule."""
    water_file = os.path.join(output_dir, "water.pdb")
    
    # Simple water molecule (O at origin, H atoms in tetrahedral positions)
    water_pdb = """HETATM    1  O   HOH     1       0.000   0.000   0.000  1.00  0.00           O
HETATM    2  H1  HOH     1       0.757   0.586   0.000  1.00  0.00           H
HETATM    3  H2  HOH     1      -0.757   0.586   0.000  1.00  0.00           H
END
"""
    
    with open(water_file, 'w') as f:
        f.write(water_pdb)
    
    return water_file


def create_packmol_input(
    protein_file: str,
    water_file: str,
    cavity_regions: List[Dict],
    output_file: str,
) -> str:
    """Create Packmol input file."""
    
    # Calculate total number of water molecules to add
    total_waters = 0
    for region in cavity_regions:
        # Estimate number of waters based on cavity volume
        n_waters = max(1, int(region["volume"] / WATERS_PER_1000_CUBIC_ANGSTROMS))
        total_waters += n_waters
    
    input_lines = [
        "tolerance 2.0",
        "filetype pdb",
        f"output {output_file}",
        "",
        "# Protein structure (fixed)",
        f"structure {protein_file}",
        "  number 1",
        "  fixed 0. 0. 0. 0. 0. 0.",
        "end structure",
        "",
    ]
    
    # Add water molecules for each cavity
    for i, region in enumerate(cavity_regions):
        n_waters = max(1, int(region["volume"] / WATERS_PER_1000_CUBIC_ANGSTROMS))
        min_coords = region["min"]
        max_coords = region["max"]
        
        input_lines.extend([
            f"# Water molecules for cavity {region['id']}",
            f"structure {water_file}",
            f"  number {n_waters}",
            f"  inside box {min_coords[0]:.3f} {min_coords[1]:.3f} {min_coords[2]:.3f} "
            f"{max_coords[0]:.3f} {max_coords[1]:.3f} {max_coords[2]:.3f}",
            "end structure",
            "",
        ])
    
    input_content = "\n".join(input_lines)
    
    # Save input file
    input_file = os.path.join(os.path.dirname(output_file), "packmol.inp")
    with open(input_file, 'w') as f:
        f.write(input_content)
    
    return input_file


def run_packmol(input_file: str, output_dir: str) -> None:
    """Run Packmol with the given input file."""
    try:
        result = subprocess.run(
            ["packmol"],
            stdin=open(input_file, 'r'),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=output_dir,
            timeout=300,  # 5 minute timeout
            text=True,
        )
        
        if result.returncode != 0:
            print(f"Packmol warning/error: {result.stderr}")
            
    except subprocess.TimeoutExpired:
        print("Warning: Packmol timed out. Using partial results if available.")
    except Exception as e:
        print(f"Warning: Error running Packmol: {e}")


def fill_with_simple_placement(
    protein_file: str,
    selected_cavities: List[Dict[str, Any]],
    cavity_data: Any,
    output_file: str,
) -> str:
    """
    Simple fallback method to place water molecules in cavities.
    Places water molecules on a grid within cavity regions.
    """
    from cavefiller.cavity_finder import get_cavity_grid_points
    
    # Read protein file
    with open(protein_file, 'r') as f:
        protein_lines = f.readlines()
    
    # Start with protein structure
    output_lines = []
    for line in protein_lines:
        if not line.startswith("END"):
            output_lines.append(line)
    
    # Add water molecules for each cavity
    water_count = 0
    for cavity in selected_cavities:
        points = get_cavity_grid_points(cavity_data, cavity["id"])
        
        if len(points) == 0:
            continue
        
        # Sample points from the cavity (not too dense)
        # Take every Nth point to avoid overcrowding
        step = max(1, len(points) // MAX_WATERS_PER_CAVITY)
        sampled_points = points[::step]
        
        for point in sampled_points:
            water_count += 1
            # Create water molecule (just oxygen for simplicity)
            x, y, z = point
            water_line = (
                f"HETATM{water_count:5d}  O   HOH  {water_count:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O\n"
            )
            output_lines.append(water_line)
    
    output_lines.append("END\n")
    
    # Write output file
    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    
    print(f"Added {water_count} water molecules using simple placement method.")
    
    return output_file
