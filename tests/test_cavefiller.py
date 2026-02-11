"""Unit tests for CaveFiller."""

import pytest
import os
import tempfile
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

EXAMPLE_PDB = Path(__file__).resolve().parents[1] / "examples" / "musM_OBP5_model_0_boltz2.pdb"


def test_imports():
    """Test that all modules can be imported."""
    from cavefiller import cavity_finder, cavity_selector, water_filler, cli
    assert cavity_finder is not None
    assert cavity_selector is not None
    assert water_filler is not None
    assert cli is not None


def test_cavity_finder_with_example_protein():
    """Test cavity finder with the Mus OBP5 example protein."""
    from cavefiller.cavity_finder import find_cavities

    if not EXAMPLE_PDB.exists():
        pytest.skip("Example protein file not found")

    with tempfile.TemporaryDirectory() as tmpdir:
        cavities, cavity_data = find_cavities(
            str(EXAMPLE_PDB),
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            output_dir=tmpdir,
        )

        # Should parse and execute without errors
        assert isinstance(cavities, list)
        # cavity_data may be None if no cavities are found


def test_read_protein_atoms():
    """Test reading protein atoms from the Mus OBP5 example PDB file."""
    from cavefiller.water_filler import read_protein_atoms

    if not EXAMPLE_PDB.exists():
        pytest.skip("Example protein file not found")

    atoms = read_protein_atoms(str(EXAMPLE_PDB))

    assert len(atoms) > 0
    # Check structure: list of (element, coordinates)
    for element, coords in atoms:
        assert isinstance(element, str)
        assert isinstance(coords, np.ndarray)
        assert len(coords) == 3


def test_clash_detection():
    """Test clash detection function."""
    from cavefiller.water_filler import check_clash
    
    # Create some protein atoms
    protein_atoms = [
        ('C', np.array([0.0, 0.0, 0.0])),
        ('N', np.array([5.0, 0.0, 0.0])),
    ]
    
    # Test position too close to protein (should clash)
    close_pos = np.array([0.5, 0.0, 0.0])
    assert check_clash(close_pos, protein_atoms, []) == True
    
    # Test position far from protein (no clash)
    far_pos = np.array([10.0, 10.0, 10.0])
    assert check_clash(far_pos, protein_atoms, []) == False
    
    # Test position close to another water (should clash)
    water_pos = np.array([3.0, 3.0, 3.0])
    existing_waters = [np.array([3.5, 3.0, 3.0])]
    assert check_clash(water_pos, [], existing_waters) == True


def test_cli_app_exists():
    """Test that the CLI app is defined."""
    from cavefiller.cli import app
    import typer
    
    assert isinstance(app, typer.Typer)


def test_cavity_with_example_protein():
    """Test cavity detection with the example protein that has a cavity."""
    from cavefiller.cavity_finder import find_cavities

    # Skip if example doesn't exist
    if not EXAMPLE_PDB.exists():
        pytest.skip("Example protein file not found")

    with tempfile.TemporaryDirectory() as tmpdir:
        cavities, cavity_data = find_cavities(
            str(EXAMPLE_PDB),
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            output_dir=tmpdir,
        )
        
        # This protein should have cavities
        assert len(cavities) > 0
        
        # Check cavity structure
        for cavity in cavities:
            assert 'id' in cavity
            assert 'volume' in cavity
            assert 'area' in cavity
            assert cavity['volume'] >= 5.0


def test_monte_carlo_water_filling():
    """Test Monte Carlo water filling method."""
    from cavefiller.water_filler import fill_cavities_with_water
    from cavefiller.cavity_finder import find_cavities

    if not EXAMPLE_PDB.exists():
        pytest.skip("Example protein file not found")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Find cavities
        cavities, cavity_data = find_cavities(
            str(EXAMPLE_PDB),
            output_dir=tmpdir,
        )
        
        if len(cavities) == 0:
            pytest.skip("No cavities found in test protein")
        
        # Fill with water using Monte Carlo
        output_file = fill_cavities_with_water(
            str(EXAMPLE_PDB),
            cavities[:1],  # Just first cavity
            cavity_data,
            tmpdir,
            waters_per_cavity={cavities[0]['id']: 5}  # Place 5 waters
        )
        
        assert os.path.exists(output_file)
        
        # Check that water was added
        with open(output_file, 'r') as f:
            content = f.read()
            assert 'HOH' in content


def test_mmff94_optimization_function_runs():
    """Test that MMFF94 fixed-protein optimization returns valid geometries."""
    from cavefiller.water_filler import optimize_waters_mmff94_fixed_protein, _protein_atoms_from_mol

    protein_mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    assert AllChem.EmbedMolecule(protein_mol, randomSeed=0xF00D) == 0
    AllChem.UFFOptimizeMolecule(protein_mol, maxIters=100)
    protein_atoms = _protein_atoms_from_mol(protein_mol)

    water_positions = [np.array([5.0, 0.0, 0.0], dtype=float)]
    water_origin_cavity_points = [
        np.array(
            [
                [4.0, 0.0, 0.0],
                [4.5, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [5.5, 0.0, 0.0],
                [6.0, 0.0, 0.0],
                [5.0, 0.5, 0.0],
                [5.0, -0.5, 0.0],
            ],
            dtype=float,
        )
    ]

    optimized_positions, optimized_geometries = optimize_waters_mmff94_fixed_protein(
        protein_mol=protein_mol,
        water_positions=water_positions,
        water_origin_cavity_points=water_origin_cavity_points,
        protein_atoms=protein_atoms,
        max_iterations=50,
    )

    assert len(optimized_positions) == 1
    assert len(optimized_geometries) == 1
    assert optimized_geometries[0][0].shape == (3,)
    assert optimized_geometries[0][1].shape == (3,)
    assert optimized_geometries[0][2].shape == (3,)


def test_fill_cavities_uses_mmff94_optimizer(monkeypatch):
    """Test that fill_cavities_with_water calls MMFF94 optimizer when enabled."""
    from cavefiller.cavity_finder import find_cavities
    from cavefiller import water_filler as wf

    if not EXAMPLE_PDB.exists():
        pytest.skip("Example protein file not found")

    with tempfile.TemporaryDirectory() as tmpdir:
        cavities, cavity_data = find_cavities(str(EXAMPLE_PDB), output_dir=tmpdir)
        if len(cavities) == 0:
            pytest.skip("No cavities found in test protein")

        called = {"value": False, "max_iterations": None}
        real_builder = wf._build_water_geometries_from_positions

        def fake_optimize(
            protein_mol,
            water_positions,
            water_origin_cavity_points,
            protein_atoms,
            max_iterations,
            remove_after_optim,
        ):
            called["value"] = True
            called["max_iterations"] = max_iterations
            return water_positions, real_builder(water_positions, protein_atoms)

        monkeypatch.setattr(wf, "optimize_waters_mmff94_fixed_protein", fake_optimize)

        output_file = wf.fill_cavities_with_water(
            str(EXAMPLE_PDB),
            cavities[:1],
            cavity_data,
            tmpdir,
            waters_per_cavity={cavities[0]["id"]: 3},
            optimize_mmff94=True,
            mmff_max_iterations=123,
        )

        assert called["value"] is True
        assert called["max_iterations"] == 123
        assert os.path.exists(output_file)


def test_fill_cavities_uses_openmm_optimizer_when_enabled(monkeypatch):
    """Test that fill_cavities_with_water calls OpenMM optimizer when enabled."""
    from cavefiller.cavity_finder import find_cavities
    from cavefiller import water_filler as wf

    if not EXAMPLE_PDB.exists():
        pytest.skip("Example protein file not found")

    with tempfile.TemporaryDirectory() as tmpdir:
        cavities, cavity_data = find_cavities(str(EXAMPLE_PDB), output_dir=tmpdir)
        if len(cavities) == 0:
            pytest.skip("No cavities found in test protein")

        called = {"value": False, "max_iterations": None}
        def fake_openmm_optimize(
            protein_mol,
            protein_file,
            water_geometries,
            protein_atoms,
            max_iterations,
            remove_after_optim,
            use_cuda,
        ):
            called["value"] = True
            called["protein_mol_atoms"] = protein_mol.GetNumAtoms()
            called["protein_file"] = protein_file
            called["max_iterations"] = max_iterations
            called["use_cuda"] = use_cuda
            return [geom[0] for geom in water_geometries], water_geometries

        def fake_mmff_optimize(*args, **kwargs):
            raise AssertionError("MMFF optimizer should not be used when OpenMM is enabled")

        monkeypatch.setattr(wf, "optimize_waters_openmm_fixed_protein", fake_openmm_optimize)
        monkeypatch.setattr(wf, "optimize_waters_mmff94_fixed_protein", fake_mmff_optimize)

        output_file = wf.fill_cavities_with_water(
            str(EXAMPLE_PDB),
            cavities[:1],
            cavity_data,
            tmpdir,
            waters_per_cavity={cavities[0]["id"]: 3},
            optimize_mmff94=True,
            optimize_openmm=True,
            openmm_max_iterations=77,
            openmm_use_cuda=True,
        )

        assert called["value"] is True
        assert called["protein_mol_atoms"] > 0
        assert called["protein_file"] == str(EXAMPLE_PDB)
        assert called["max_iterations"] == 77
        assert called["use_cuda"] is True
        assert os.path.exists(output_file)


def test_fill_cavities_keep_all_sets_mmff_drop_behavior(monkeypatch):
    """Test that keep_all toggles MMFF post-optimization dropping behavior."""
    from cavefiller.cavity_finder import find_cavities
    from cavefiller import water_filler as wf

    if not EXAMPLE_PDB.exists():
        pytest.skip("Example protein file not found")

    with tempfile.TemporaryDirectory() as tmpdir:
        cavities, cavity_data = find_cavities(str(EXAMPLE_PDB), output_dir=tmpdir)
        if len(cavities) == 0:
            pytest.skip("No cavities found in test protein")

        called = {"remove_after_optim": None}
        real_builder = wf._build_water_geometries_from_positions

        def fake_optimize(
            protein_mol,
            water_positions,
            water_origin_cavity_points,
            protein_atoms,
            max_iterations,
            remove_after_optim,
        ):
            called["remove_after_optim"] = remove_after_optim
            return water_positions, real_builder(water_positions, protein_atoms)

        monkeypatch.setattr(wf, "optimize_waters_mmff94_fixed_protein", fake_optimize)

        output_file = wf.fill_cavities_with_water(
            str(EXAMPLE_PDB),
            cavities[:1],
            cavity_data,
            tmpdir,
            waters_per_cavity={cavities[0]["id"]: 2},
            optimize_mmff94=True,
            keep_all=True,
        )

        assert called["remove_after_optim"] is False


def test_local_relaxation_separates_overlapping_waters():
    """Local relaxation should push colliding waters apart."""
    from cavefiller import water_filler as wf

    protein_atoms = [("C", np.array([0.0, 3.2, 0.0], dtype=float))]
    initial_positions = [
        np.array([0.0, 0.0, 0.0], dtype=float),
        np.array([0.2, 0.0, 0.0], dtype=float),
    ]
    geometries = wf._build_water_geometries_from_positions(initial_positions, protein_atoms)

    cavity_track = np.array([[x, 0.0, 0.0] for x in np.linspace(-3.0, 3.0, 61)], dtype=float)
    origins = [cavity_track, cavity_track]

    relaxed = wf._relax_water_geometries_locally(geometries, origins, protein_atoms)
    relaxed_o0 = relaxed[0][0]
    relaxed_o1 = relaxed[1][0]
    relaxed_dist = np.linalg.norm(relaxed_o0 - relaxed_o1)

    assert relaxed_dist > np.linalg.norm(initial_positions[0] - initial_positions[1])
    assert relaxed_dist >= 2.0


def test_mmff_selection_prefers_higher_clearance_candidate():
    """Selection should choose valid candidate with better all-atom clearance."""
    from cavefiller import water_filler as wf

    protein_atoms = [("C", np.array([0.0, 0.0, 0.0], dtype=float))]
    optimized_geom = wf._build_water_geometries_from_positions(
        [np.array([2.8, 0.0, 0.0], dtype=float)],
        protein_atoms,
    )[0]
    original_geom = wf._build_water_geometries_from_positions(
        [np.array([4.0, 0.0, 0.0], dtype=float)],
        protein_atoms,
    )[0]
    origin_cavity = [np.array([[2.8, 0.0, 0.0], [4.0, 0.0, 0.0]], dtype=float)]

    selected = wf._select_valid_geometries_after_mmff(
        optimized_geometries=[optimized_geom],
        original_geometries=[original_geom],
        water_origin_cavity_points=origin_cavity,
        protein_atoms=protein_atoms,
    )

    assert len(selected) == 1
    assert np.allclose(selected[0][0], original_geom[0])
