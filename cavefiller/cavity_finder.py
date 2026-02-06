"""Cavity detection using pyKVFinder."""

import os
from typing import List, Dict, Tuple, Any
import numpy as np

# Grid spacing for cavity detection (in Angstroms)
DEFAULT_GRID_STEP = 0.6


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
        volume_cutoff: Minimum cavity volume to consider (Ų)
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
        input=protein_file,
        probe_in=probe_in,
        probe_out=probe_out,
        step=DEFAULT_GRID_STEP,  # Grid step size
        volume_cutoff=volume_cutoff,
    )
    
    # Extract cavity information
    cavities = []
    
    # Get cavity volumes and areas
    if hasattr(cavity_data, 'volume') and cavity_data.volume is not None:
        volumes = cavity_data.volume
        areas = cavity_data.area if hasattr(cavity_data, 'area') else {}
        
        # Create mapping from string IDs to integer IDs
        # KVFinder uses string IDs like 'KAA', 'KAB', etc., but the grid uses integers
        cavity_id_map = {}
        for idx, (cavity_str_id, volume) in enumerate(volumes.items(), start=1):
            cavity_id_map[cavity_str_id] = idx
        
        # Process each cavity
        for cavity_str_id, volume in volumes.items():
            if volume >= volume_cutoff:
                cavity_info = {
                    "id": cavity_id_map[cavity_str_id],
                    "string_id": cavity_str_id,
                    "volume": volume,
                    "area": areas.get(cavity_str_id, 0.0) if areas else 0.0,
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
        cavity_id: Integer ID of the cavity (1-indexed)
        
    Returns:
        Array of (x, y, z) coordinates for the cavity grid points
    """
    if not hasattr(cavity_data, 'cavities') or cavity_data.cavities is None:
        return np.array([])
    
    # Get cavity grid
    cavity_grid = cavity_data.cavities
    
    # Find all points belonging to this cavity
    # Note: KVFinder uses 1-indexed cavity IDs in the grid
    points = np.argwhere(cavity_grid == cavity_id)
    
    # Convert grid indices to real coordinates if origin metadata is available.
    # Different pyKVFinder versions expose metadata on either cavity_data or cavity_data.surface.
    step = getattr(cavity_data, "step", DEFAULT_GRID_STEP)
    origin = None

    if hasattr(cavity_data, "surface") and hasattr(cavity_data.surface, "P1"):
        origin = np.array([cavity_data.surface.P1[i] for i in range(3)], dtype=float)
    elif hasattr(cavity_data, "P1"):
        origin = np.array([cavity_data.P1[i] for i in range(3)], dtype=float)

    points = points.astype(float)
    if origin is not None:
        return origin + points * float(step)

    # Fallback: return index-space points; downstream code will align to protein frame.
    return points
