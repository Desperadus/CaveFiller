"""Cavity detection using pyKVFinder."""

import os
from typing import List, Dict, Tuple, Any
import numpy as np


def find_cavities(
    protein_file: str,
    probe_in: float = 1.4,
    probe_out: float = 4.0,
    volume_cutoff: float = 5.0,
    output_dir: str = "./output",
) -> Tuple[List[Dict[str, Any]], Any]:
    """
    Find cavities in a protein structure using pyKVFinder.
    
    Args:
        protein_file: Path to the protein PDB file
        probe_in: Probe In radius for cavity detection (Å)
        probe_out: Probe Out radius for cavity detection (Å)
        volume_cutoff: Minimum cavity volume to consider (Å³)
        output_dir: Directory to save cavity detection results
        
    Returns:
        Tuple of (list of cavity dictionaries, cavity_data object)
    """
    try:
        import pyKVFinder
    except ImportError:
        raise ImportError(
            "pyKVFinder is not installed. Please install it with: pip install pykvfinder"
        )
    
    # Run KVFinder to detect cavities
    cavity_data = pyKVFinder.run_workflow(
        pdb=protein_file,
        probe_in=probe_in,
        probe_out=probe_out,
        step=0.6,  # Grid step size
        output=os.path.join(output_dir, "cavities.toml"),
    )
    
    # Extract cavity information
    cavities = []
    
    # Get cavity volumes and areas
    if hasattr(cavity_data, 'volume') and cavity_data.volume is not None:
        volumes = cavity_data.volume
        areas = cavity_data.area if hasattr(cavity_data, 'area') else {}
        
        # Process each cavity
        for cavity_id, volume in volumes.items():
            if volume >= volume_cutoff:
                cavity_info = {
                    "id": cavity_id,
                    "volume": volume,
                    "area": areas.get(cavity_id, 0.0) if areas else 0.0,
                }
                cavities.append(cavity_info)
    
    # Sort cavities by volume (largest first)
    cavities.sort(key=lambda x: x["volume"], reverse=True)
    
    return cavities, cavity_data


def get_cavity_grid_points(cavity_data: Any, cavity_id: int) -> np.ndarray:
    """
    Get the grid points that belong to a specific cavity.
    
    Args:
        cavity_data: The cavity data object from pyKVFinder
        cavity_id: ID of the cavity
        
    Returns:
        Array of (x, y, z) coordinates for the cavity grid points
    """
    if not hasattr(cavity_data, 'cavities') or cavity_data.cavities is None:
        return np.array([])
    
    # Get cavity grid
    cavity_grid = cavity_data.cavities
    
    # Find all points belonging to this cavity
    points = np.argwhere(cavity_grid == cavity_id)
    
    # Convert grid indices to real coordinates
    if hasattr(cavity_data, 'origin') and hasattr(cavity_data, 'step'):
        origin = np.array(cavity_data.origin)
        step = cavity_data.step
        real_coords = origin + points * step
        return real_coords
    
    return points
