"""Fill selected cavities with explicit water molecules using RDKit."""

import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from cavefiller.cavity_finder import get_cavity_grid_points

# Van der Waals radii for clash detection (Angstrom)
VDW_RADII = {
    "H": 1.2,
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "S": 1.8,
    "P": 1.8,
}

# Conservative spacing to avoid overlaps
MIN_WATER_WATER_DISTANCE = 2.7
MIN_WATER_PROTEIN_DISTANCE = 2.35
MAX_WATER_PROTEIN_DISTANCE = 5.5

# Water geometry
OH_BOND_LENGTH = 0.9572
HOH_ANGLE_DEG = 104.52
DEFAULT_MMFF_MAX_ITERS = 300


def fill_cavities_with_water(
    protein_file: str,
    selected_cavities: List[Dict[str, Any]],
    cavity_data: Any,
    output_dir: str,
    waters_per_cavity: Dict[int, int] = None,
    optimize_mmff94: bool = True,
    mmff_max_iterations: int = DEFAULT_MMFF_MAX_ITERS,
) -> str:
    """Place explicit waters in selected cavities and write a combined PDB."""
    base_name = os.path.splitext(os.path.basename(protein_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}_filled.pdb")
    optimized_output_file = os.path.join(output_dir, f"{base_name}_filled_optimized.pdb")
    protein_mol = _load_protein_mol(protein_file)
    protein_atoms = read_protein_atoms(protein_file)

    all_water_positions: List[np.ndarray] = []
    all_water_origin_cavity_points: List[np.ndarray] = []

    for cavity in selected_cavities:
        cavity_id = cavity["id"]
        volume = cavity["volume"]

        if waters_per_cavity and cavity_id in waters_per_cavity:
            n_waters = waters_per_cavity[cavity_id]
        else:
            n_waters = max(1, int(volume / 30.0))

        print(
            f"Placing up to {n_waters} waters in cavity {cavity_id} "
            f"(volume: {volume:.2f} A^3)..."
        )

        cavity_points = get_cavity_grid_points(cavity_data, cavity_id)
        if len(cavity_points) == 0:
            print(f"  Warning: no grid points found for cavity {cavity_id}")
            continue

        cavity_points = _ensure_cavity_points_near_protein(cavity_points, protein_atoms)

        positions = monte_carlo_water_placement(
            cavity_points,
            protein_atoms,
            n_waters,
            max_attempts=500,
        )

        if positions:
            all_water_positions.extend(positions)
            all_water_origin_cavity_points.extend([cavity_points] * len(positions))
            print(f"  Placed {len(positions)} waters")
        else:
            print("  Warning: no waters could be placed")

    if all_water_positions:
        nonoptimized_geometries = _build_water_geometries_from_positions(
            all_water_positions,
            protein_atoms,
        )
        waters_mol = build_waters_mol(
            water_positions=[],
            chain_id="W",
            water_geometries=nonoptimized_geometries,
        )
        combined = Chem.CombineMols(protein_mol, waters_mol)
        Chem.MolToPDBFile(combined, output_file)
        print(f"Saved non-optimized structure: {output_file}")

        if optimize_mmff94 and all_water_origin_cavity_points:
            optimized_positions, optimized_geometries = optimize_waters_mmff94_fixed_protein(
                protein_mol=protein_mol,
                water_positions=all_water_positions,
                water_origin_cavity_points=all_water_origin_cavity_points,
                protein_atoms=protein_atoms,
                max_iterations=mmff_max_iterations,
            )
            if optimized_positions:
                print(
                    f"MMFF94 (fixed protein) kept {len(optimized_positions)}/"
                    f"{len(all_water_positions)} waters after filtering"
                )
                optimized_waters_mol = build_waters_mol(
                    water_positions=[],
                    chain_id="W",
                    water_geometries=optimized_geometries,
                )
                optimized_combined = Chem.CombineMols(protein_mol, optimized_waters_mol)
                Chem.MolToPDBFile(optimized_combined, optimized_output_file)
                print(f"Saved optimized structure: {optimized_output_file}")
                all_water_positions = optimized_positions
                output_file = optimized_output_file
            else:
                print("Warning: MMFF94 step yielded no valid waters, keeping initial placement")
    else:
        Chem.MolToPDBFile(protein_mol, output_file)

    print(f"\nTotal waters added: {len(all_water_positions)}")
    return output_file


def _load_protein_mol(protein_file: str) -> Chem.Mol:
    """Load protein PDB into RDKit while preserving atom records."""
    mol = Chem.MolFromPDBFile(protein_file, sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not read protein file as PDB: {protein_file}")
    return mol


def read_protein_atoms(protein_file: str) -> List[Tuple[str, np.ndarray]]:
    """Read protein atom element and coordinates from RDKit PDB parsing."""
    mol = _load_protein_mol(protein_file)
    conf = mol.GetConformer()

    atoms: List[Tuple[str, np.ndarray]] = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms.append(
            (
                atom.GetSymbol(),
                np.array([float(pos.x), float(pos.y), float(pos.z)], dtype=float),
            )
        )
    return atoms


def check_clash(
    position: np.ndarray,
    existing_atoms: List[Tuple[str, np.ndarray]],
    existing_waters: List[np.ndarray],
) -> bool:
    """Return True when a candidate oxygen clashes with protein or waters."""
    for element, atom_coords in existing_atoms:
        distance = np.linalg.norm(position - atom_coords)
        vdw_radius = VDW_RADII.get(element, 1.7)
        min_distance = max(
            vdw_radius + VDW_RADII["O"] - 0.5,
            MIN_WATER_PROTEIN_DISTANCE,
        )
        if distance < min_distance:
            return True

    for water_pos in existing_waters:
        if np.linalg.norm(position - water_pos) < MIN_WATER_WATER_DISTANCE:
            return True

    return False


def monte_carlo_water_placement(
    cavity_points: np.ndarray,
    protein_atoms: List[Tuple[str, np.ndarray]],
    n_waters: int,
    max_attempts: int = 500,
) -> List[np.ndarray]:
    """Place waters by sampling cavity grid points plus small local jitter."""
    if len(cavity_points) == 0 or n_waters <= 0:
        return []

    placed_waters: List[np.ndarray] = []
    cavity_points = np.asarray(cavity_points, dtype=float)
    protein_coords = np.asarray([coords for _, coords in protein_atoms], dtype=float)

    for _ in range(n_waters):
        placed = False
        for _attempt in range(max_attempts):
            base = cavity_points[np.random.randint(0, len(cavity_points))]
            jitter = np.random.normal(loc=0.0, scale=0.22, size=3)
            position = base + jitter

            # Keep candidate in/near cavity voxels only.
            nearest = np.min(np.linalg.norm(cavity_points - position, axis=1))
            if nearest > 0.7:
                continue

            # Waters must remain physically close to the protein envelope.
            nearest_protein = np.min(np.linalg.norm(protein_coords - position, axis=1))
            if nearest_protein > MAX_WATER_PROTEIN_DISTANCE:
                continue

            if not check_clash(position, protein_atoms, placed_waters):
                placed_waters.append(position)
                placed = True
                break

        if not placed:
            print(f"    Warning: placed {len(placed_waters)}/{n_waters} waters")
            break

    return placed_waters


def _protein_atoms_from_mol(mol: Chem.Mol) -> List[Tuple[str, np.ndarray]]:
    """Extract (element, coords) from an RDKit molecule conformer."""
    conf = mol.GetConformer()
    atoms: List[Tuple[str, np.ndarray]] = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms.append(
            (
                atom.GetSymbol(),
                np.array([float(pos.x), float(pos.y), float(pos.z)], dtype=float),
            )
        )
    return atoms


def _prepare_protein_for_mmff(protein_mol: Chem.Mol) -> Chem.Mol:
    """Prepare protein for MMFF and attempt explicit hydrogenation."""
    mol = Chem.Mol(protein_mol)
    mol.UpdatePropertyCache(strict=False)
    try:
        Chem.GetSymmSSSR(mol)
    except Exception:
        pass
    try:
        AllChem.MMFFSanitizeMolecule(mol)
    except Exception:
        pass
    try:
        mol = Chem.AddHs(mol, addCoords=True)
    except Exception as exc:
        print(f"Warning: could not add explicit protein hydrogens for MMFF ({exc})")
    return mol


def _principal_direction_from_protein(
    oxygen: np.ndarray,
    protein_coords: Optional[np.ndarray],
) -> np.ndarray:
    """Direction from nearest protein atom to oxygen (points away from protein)."""
    if protein_coords is None or len(protein_coords) == 0:
        return np.array([1.0, 0.0, 0.0], dtype=float)
    deltas = oxygen - protein_coords
    dists = np.linalg.norm(deltas, axis=1)
    nearest = int(np.argmin(dists))
    vec = deltas[nearest]
    norm = np.linalg.norm(vec)
    if norm < 1e-8:
        return np.array([1.0, 0.0, 0.0], dtype=float)
    return vec / norm


def optimize_waters_mmff94_fixed_protein(
    protein_mol: Chem.Mol,
    water_positions: List[np.ndarray],
    water_origin_cavity_points: List[np.ndarray],
    protein_atoms: List[Tuple[str, np.ndarray]],
    max_iterations: int = DEFAULT_MMFF_MAX_ITERS,
) -> Tuple[List[np.ndarray], List[Tuple[np.ndarray, np.ndarray, np.ndarray]]]:
    """MMFF94 optimize waters with protein atoms fixed in place."""
    if not water_positions:
        return [], []

    try:
        protein_mol_for_mmff = _prepare_protein_for_mmff(protein_mol)
        protein_atoms_for_mmff = _protein_atoms_from_mol(protein_mol_for_mmff)
        original_geometries = _build_water_geometries_from_positions(
            water_positions, protein_atoms_for_mmff
        )
        waters_mol = build_waters_mol(
            water_positions=[],
            chain_id="W",
            water_geometries=original_geometries,
        )
        complex_mol = Chem.CombineMols(protein_mol_for_mmff, waters_mol)
        work_mol = Chem.Mol(complex_mol)

        # MMFF expects initialized ring/aromaticity/cache data.
        work_mol.UpdatePropertyCache(strict=False)
        try:
            Chem.GetSymmSSSR(work_mol)
        except Exception:
            pass
        try:
            AllChem.MMFFSanitizeMolecule(work_mol)
        except Exception:
            # Continue and let MMFF parameter checks decide viability.
            pass

        if not AllChem.MMFFHasAllMoleculeParams(work_mol):
            print("Warning: MMFF94 params missing for complex; skipping MMFF optimization")
            return water_positions, original_geometries

        props = AllChem.MMFFGetMoleculeProperties(work_mol, mmffVariant="MMFF94")
        if props is None:
            print("Warning: MMFF94 property setup failed; skipping MMFF optimization")
            return water_positions, original_geometries

        ff = AllChem.MMFFGetMoleculeForceField(work_mol, props, confId=0)
        if ff is None:
            print("Warning: MMFF94 force field setup failed; skipping MMFF optimization")
            return water_positions, original_geometries

        protein_atom_count = protein_mol_for_mmff.GetNumAtoms()
        for atom_idx in range(protein_atom_count):
            ff.AddFixedPoint(atom_idx)

        ff.Initialize()
        ff.Minimize(maxIts=max_iterations)

        conf = work_mol.GetConformer()
        optimized_geometries = []
        for i in range(len(water_positions)):
            oxygen_idx = protein_atom_count + i * 3
            h1_idx = oxygen_idx + 1
            h2_idx = oxygen_idx + 2
            o = conf.GetAtomPosition(oxygen_idx)
            h1 = conf.GetAtomPosition(h1_idx)
            h2 = conf.GetAtomPosition(h2_idx)
            optimized_geometries.append(
                (
                    np.array([float(o.x), float(o.y), float(o.z)], dtype=float),
                    np.array([float(h1.x), float(h1.y), float(h1.z)], dtype=float),
                    np.array([float(h2.x), float(h2.y), float(h2.z)], dtype=float),
                )
            )

        selected_geometries = _select_valid_geometries_after_mmff(
            optimized_geometries=optimized_geometries,
            original_geometries=original_geometries,
            water_origin_cavity_points=water_origin_cavity_points,
            protein_atoms=protein_atoms_for_mmff,
        )
        return [geom[0] for geom in selected_geometries], selected_geometries
    except Exception as exc:
        print(f"Warning: MMFF94 optimization failed ({exc}); keeping initial placement")
        fallback_geometries = _build_water_geometries_from_positions(
            water_positions, protein_atoms
        )
        return water_positions, fallback_geometries


def _is_position_valid_for_origin(
    position: np.ndarray,
    cavity_points: np.ndarray,
    protein_atoms: List[Tuple[str, np.ndarray]],
) -> bool:
    """Check if a water oxygen is valid for its original cavity."""
    cavity_points = np.asarray(cavity_points, dtype=float)
    protein_coords = np.asarray([coords for _, coords in protein_atoms], dtype=float)

    nearest_cavity = np.min(np.linalg.norm(cavity_points - position, axis=1))
    if nearest_cavity > 1.0:
        return False

    nearest_protein = np.min(np.linalg.norm(protein_coords - position, axis=1))
    if nearest_protein > MAX_WATER_PROTEIN_DISTANCE:
        return False

    if check_clash(position, protein_atoms, []):
        return False

    return True


def _min_allowed_pair_distance(elem_a: str, elem_b: str) -> float:
    """Conservative lower-bound distance for two atoms."""
    vdw_a = VDW_RADII.get(elem_a, 1.7)
    vdw_b = VDW_RADII.get(elem_b, 1.7)
    # Allow tighter contacts when hydrogen is involved, but prevent hard overlaps.
    tol = 1.0 if "H" in (elem_a, elem_b) else 0.5
    raw = vdw_a + vdw_b - tol
    if elem_a == "H" and elem_b == "H":
        return max(1.35, raw)
    if "H" in (elem_a, elem_b):
        return max(1.5, raw)
    return max(2.0, raw)


def _geometry_has_all_atom_clash(
    geometry: Tuple[np.ndarray, np.ndarray, np.ndarray],
    protein_atoms: List[Tuple[str, np.ndarray]],
    accepted_geometries: List[Tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> bool:
    """All-atom clash check for one water geometry."""
    water_atoms = [
        ("O", geometry[0]),
        ("H", geometry[1]),
        ("H", geometry[2]),
    ]

    for w_elem, w_pos in water_atoms:
        for p_elem, p_pos in protein_atoms:
            if np.linalg.norm(w_pos - p_pos) < _min_allowed_pair_distance(w_elem, p_elem):
                return True

    for accepted in accepted_geometries:
        accepted_atoms = [
            ("O", accepted[0]),
            ("H", accepted[1]),
            ("H", accepted[2]),
        ]
        for w_elem, w_pos in water_atoms:
            for a_elem, a_pos in accepted_atoms:
                if np.linalg.norm(w_pos - a_pos) < _min_allowed_pair_distance(w_elem, a_elem):
                    return True

    return False


def _select_valid_geometries_after_mmff(
    optimized_geometries: List[Tuple[np.ndarray, np.ndarray, np.ndarray]],
    original_geometries: List[Tuple[np.ndarray, np.ndarray, np.ndarray]],
    water_origin_cavity_points: List[np.ndarray],
    protein_atoms: List[Tuple[str, np.ndarray]],
) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Keep each water in original cavity after MMFF; enforce all-atom clash checks."""
    accepted: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    reverted_count = 0
    dropped_count = 0

    for optimized_geom, original_geom, origin_cavity in zip(
        optimized_geometries, original_geometries, water_origin_cavity_points
    ):
        candidates: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []

        if _is_position_valid_for_origin(optimized_geom[0], origin_cavity, protein_atoms):
            candidates.append(optimized_geom)
        else:
            reverted_count += 1

        if _is_position_valid_for_origin(original_geom[0], origin_cavity, protein_atoms):
            candidates.append(original_geom)

        chosen = None
        for candidate in candidates:
            if not _geometry_has_all_atom_clash(candidate, protein_atoms, accepted):
                chosen = candidate
                break

        if chosen is None:
            dropped_count += 1
            continue

        accepted.append(chosen)

    if reverted_count:
        print(
            f"MMFF94 check: reverted {reverted_count} waters that moved out "
            "of their original cavity"
        )
    if dropped_count:
        print(f"MMFF94 check: dropped {dropped_count} waters after validation/clash checks")

    return accepted


def _points_overlap_ratio(points: np.ndarray, box_min: np.ndarray, box_max: np.ndarray) -> float:
    """Fraction of points inside a given axis-aligned box."""
    inside = np.all((points >= box_min) & (points <= box_max), axis=1)
    return float(np.mean(inside)) if len(points) else 0.0


def _ensure_cavity_points_near_protein(
    cavity_points: np.ndarray, protein_atoms: List[Tuple[str, np.ndarray]]
) -> np.ndarray:
    """
    Align cavity points to protein frame when KVFinder points look off-frame.

    Some KVFinder runs may expose cavity grids in a different axis/origin frame.
    This picks the best of common transforms by maximizing overlap with protein bbox.
    """
    points = np.asarray(cavity_points, dtype=float)
    if len(points) == 0:
        return points

    protein_coords = np.asarray([coords for _, coords in protein_atoms], dtype=float)
    pmin = protein_coords.min(axis=0) - 5.0
    pmax = protein_coords.max(axis=0) + 5.0
    protein_centroid = protein_coords.mean(axis=0)

    base_candidates = [points, points[:, [2, 1, 0]]]
    candidates = []
    for base in base_candidates:
        candidates.append(base)
        centroid_shift = protein_centroid - base.mean(axis=0)
        candidates.append(base + centroid_shift)

    best = points
    best_score = _points_overlap_ratio(points, pmin, pmax)
    for cand in candidates:
        score = _points_overlap_ratio(cand, pmin, pmax)
        if score > best_score:
            best_score = score
            best = cand

    if best_score < 0.05:
        print(
            "  Warning: cavity grid appears poorly aligned with protein coordinates; "
            "water placement may be limited"
        )

    return best


def _orthonormal_perpendicular(u: np.ndarray) -> np.ndarray:
    """Return a deterministic unit vector perpendicular to u."""
    if abs(float(u[0])) < 0.9:
        ref = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)
    v = np.cross(u, ref)
    n = np.linalg.norm(v)
    if n < 1e-8:
        ref = np.array([0.0, 0.0, 1.0], dtype=float)
        v = np.cross(u, ref)
        n = np.linalg.norm(v)
    if n < 1e-8:
        return np.array([0.0, 1.0, 0.0], dtype=float)
    return v / n


def _build_water_geometries_from_positions(
    water_positions: List[np.ndarray],
    protein_atoms: List[Tuple[str, np.ndarray]],
) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Build deterministic H-O-H geometries from oxygen positions."""
    angle_rad = np.deg2rad(HOH_ANGLE_DEG)
    cos_t = np.cos(angle_rad)
    sin_t = np.sin(angle_rad)
    protein_coords = np.asarray([coords for _, coords in protein_atoms], dtype=float)

    geometries: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for oxygen in water_positions:
        o = np.asarray(oxygen, dtype=float)
        u = _principal_direction_from_protein(o, protein_coords)
        v = _orthonormal_perpendicular(u)
        h1 = o + OH_BOND_LENGTH * u
        h2 = o + OH_BOND_LENGTH * (cos_t * u + sin_t * v)
        geometries.append((o, h1, h2))
    return geometries


def build_waters_mol(
    water_positions: List[np.ndarray],
    chain_id: str = "W",
    water_geometries: Optional[List[Tuple[np.ndarray, np.ndarray, np.ndarray]]] = None,
) -> Chem.Mol:
    """Create an RDKit molecule containing explicit HOH residues."""
    if water_geometries is None:
        raise ValueError("water_geometries must be provided")
    if not water_geometries:
        return Chem.Mol()

    water_template = Chem.AddHs(Chem.MolFromSmiles("O"))
    waters_mol = None

    for idx, geometry in enumerate(water_geometries, start=1):
        oxygen, h1, h2 = geometry

        mol = Chem.Mol(water_template)
        conf = Chem.Conformer(mol.GetNumAtoms())
        conf.SetAtomPosition(0, (float(oxygen[0]), float(oxygen[1]), float(oxygen[2])))
        conf.SetAtomPosition(1, (float(h1[0]), float(h1[1]), float(h1[2])))
        conf.SetAtomPosition(2, (float(h2[0]), float(h2[1]), float(h2[2])))
        mol.RemoveAllConformers()
        mol.AddConformer(conf)

        atom_labels = [" O  ", " H1 ", " H2 "]
        for atom_i, atom in enumerate(mol.GetAtoms()):
            info = Chem.AtomPDBResidueInfo()
            info.SetName(atom_labels[atom_i])
            info.SetResidueName("HOH")
            info.SetResidueNumber(idx)
            info.SetChainId(chain_id)
            info.SetIsHeteroAtom(True)
            atom.SetMonomerInfo(info)

        waters_mol = mol if waters_mol is None else Chem.CombineMols(waters_mol, mol)

    return waters_mol
