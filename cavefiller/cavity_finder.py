"""Cavity detection using pyKVFinder."""

from typing import List, Dict, Tuple, Any
import numpy as np

# Grid spacing for cavity detection (in Angstroms)
DEFAULT_GRID_STEP = 0.6
DEFAULT_PROBE_IN = 1.4
DEFAULT_PROBE_OUT = 4.0
DEFAULT_EXTERIOR_TRIM_DISTANCE = 2.4
DEFAULT_VOLUME_CUTOFF = 5.0


def _map_volume_keys_to_grid_labels(cavity_data: Any, volume_keys: List[str]) -> Dict[str, int]:
    """
    Map pyKVFinder cavity string IDs (KAA, KAB, ...) to integer labels in `cavity_data.cavities`.

    pyKVFinder 0.9.0 often uses positive cavity labels starting at 2 (label 1 is reserved),
    while volumes are reported by sequential string IDs. We map by ordered positive labels.
    """
    labels = sorted(int(v) for v in np.unique(cavity_data.cavities) if int(v) > 0)
    if not labels:
        return {}

    # Most pyKVFinder outputs are contiguous and aligned with volume-key order.
    if len(labels) >= len(volume_keys):
        ordered_labels = labels[: len(volume_keys)]
    else:
        ordered_labels = labels + list(range(labels[-1] + 1, labels[-1] + 1 + (len(volume_keys) - len(labels))))

    return {key: int(ordered_labels[idx]) for idx, key in enumerate(volume_keys)}


def find_cavities(
    protein_file: str,
    probe_in: float = DEFAULT_PROBE_IN,
    probe_out: float = DEFAULT_PROBE_OUT,
    step: float = DEFAULT_GRID_STEP,
    removal_distance: float = DEFAULT_EXTERIOR_TRIM_DISTANCE,
    volume_cutoff: float = DEFAULT_VOLUME_CUTOFF,
    output_dir: str = "./output",
) -> Tuple[List[Dict[str, Any]], Any]:
    """
    Find cavities in a protein structure using pyKVFinder.
    
    Args:
        protein_file: Path to the protein PDB file
        probe_in: Probe In radius for cavity detection (Å)
        probe_out: Probe Out radius for cavity detection (Å)
        step: Grid spacing for cavity detection (Å)
        removal_distance: Exterior trim distance for cavity detection (Å)
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
        step=step,
        removal_distance=removal_distance,
        volume_cutoff=volume_cutoff,
    )
    
    # Extract cavity information
    cavities = []
    
    # Get cavity volumes and areas
    if hasattr(cavity_data, 'volume') and cavity_data.volume is not None:
        volumes = cavity_data.volume
        areas = cavity_data.area if hasattr(cavity_data, 'area') else {}
        
        # Map string IDs to underlying integer cavity-grid labels.
        # User-facing IDs remain sequential for compatibility.
        volume_keys = list(volumes.keys())
        cavity_grid_id_map = _map_volume_keys_to_grid_labels(cavity_data, volume_keys)

        # Process each cavity
        for display_idx, (cavity_str_id, volume) in enumerate(volumes.items(), start=1):
            if volume >= volume_cutoff:
                cavity_info = {
                    "id": display_idx,
                    "grid_id": cavity_grid_id_map.get(cavity_str_id, display_idx),
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
    
    points = points.astype(float)

    # Convert grid indices to real coordinates if metadata is available.
    # pyKVFinder versions expose this either as public P1/P2/P3/P4 or private _vertices/_step.
    step = float(getattr(cavity_data, "step", getattr(cavity_data, "_step", DEFAULT_GRID_STEP)))
    vertices = None

    if hasattr(cavity_data, "surface") and hasattr(cavity_data.surface, "P1"):
        vertices = np.array(
            [
                [cavity_data.surface.P1[i] for i in range(3)],
                [cavity_data.surface.P2[i] for i in range(3)],
                [cavity_data.surface.P3[i] for i in range(3)],
                [cavity_data.surface.P4[i] for i in range(3)],
            ],
            dtype=float,
        )
    elif hasattr(cavity_data, "P1"):
        vertices = np.array(
            [
                [cavity_data.P1[i] for i in range(3)],
                [cavity_data.P2[i] for i in range(3)],
                [cavity_data.P3[i] for i in range(3)],
                [cavity_data.P4[i] for i in range(3)],
            ],
            dtype=float,
        )
    elif hasattr(cavity_data, "_vertices"):
        vertices = np.asarray(cavity_data._vertices, dtype=float)

    if vertices is not None and vertices.shape[0] >= 4:
        origin = vertices[0]
        axes = [vertices[1] - origin, vertices[2] - origin, vertices[3] - origin]
        unit_axes = []
        for axis in axes:
            norm = np.linalg.norm(axis)
            unit_axes.append(axis / norm if norm > 1e-8 else np.zeros(3, dtype=float))
        return (
            origin
            + points[:, [0]] * (step * unit_axes[0])
            + points[:, [1]] * (step * unit_axes[1])
            + points[:, [2]] * (step * unit_axes[2])
        )
    if vertices is not None and vertices.shape[0] >= 1:
        return vertices[0] + points * step

    # Fallback: return index-space points; downstream code will align to protein frame.
    return points
