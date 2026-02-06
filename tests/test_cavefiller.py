"""Unit tests for CaveFiller."""

import pytest
import os
import tempfile
from pathlib import Path
import numpy as np


def create_simple_protein_pdb(filepath):
    """Create a simple protein PDB file for testing."""
    content = """HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.258   2.390   0.000  1.00  0.00           O
END
"""
    with open(filepath, 'w') as f:
        f.write(content)


def test_imports():
    """Test that all modules can be imported."""
    from cavefiller import cavity_finder, cavity_selector, water_filler, cli
    assert cavity_finder is not None
    assert cavity_selector is not None
    assert water_filler is not None
    assert cli is not None


def test_cavity_finder_with_small_protein():
    """Test cavity finder with a small protein (may not find cavities)."""
    from cavefiller.cavity_finder import find_cavities
    
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_file = os.path.join(tmpdir, "test.pdb")
        create_simple_protein_pdb(pdb_file)
        
        cavities, cavity_data = find_cavities(
            pdb_file,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            output_dir=tmpdir,
        )
        
        # Small protein may not have cavities, but should not error
        assert isinstance(cavities, list)
        # cavity_data may be None if no cavities are found


def test_read_protein_atoms():
    """Test reading protein atoms from PDB file."""
    from cavefiller.water_filler import read_protein_atoms
    
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_file = os.path.join(tmpdir, "test.pdb")
        create_simple_protein_pdb(pdb_file)
        
        atoms = read_protein_atoms(pdb_file)
        
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


def test_monte_carlo_placement():
    """Test Monte Carlo water placement."""
    from cavefiller.water_filler import monte_carlo_water_placement
    
    # Create a simple cavity (grid of points)
    cavity_points = np.array([
        [10.0, 10.0, 10.0],
        [10.5, 10.0, 10.0],
        [10.0, 10.5, 10.0],
        [10.0, 10.0, 10.5],
        [11.0, 10.0, 10.0],
        [10.0, 11.0, 10.0],
        [10.0, 10.0, 11.0],
    ])
    
    # No protein atoms nearby
    protein_atoms = [
        ('C', np.array([0.0, 0.0, 0.0])),
    ]
    
    # Try to place 2 waters
    waters = monte_carlo_water_placement(cavity_points, protein_atoms, 2, max_attempts=100)
    
    # Should be able to place at least one water
    assert len(waters) >= 1
    assert len(waters) <= 2
    
    # Check that waters are numpy arrays
    for water in waters:
        assert isinstance(water, np.ndarray)
        assert len(water) == 3


def test_cli_app_exists():
    """Test that the CLI app is defined."""
    from cavefiller.cli import app
    import typer
    
    assert isinstance(app, typer.Typer)


def test_cavity_with_example_protein():
    """Test cavity detection with the example protein that has a cavity."""
    from cavefiller.cavity_finder import find_cavities
    
    example_pdb = "examples/protein_with_cavity.pdb"
    
    # Skip if example doesn't exist
    if not os.path.exists(example_pdb):
        pytest.skip("Example protein file not found")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        cavities, cavity_data = find_cavities(
            example_pdb,
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
    
    example_pdb = "examples/protein_with_cavity.pdb"
    
    if not os.path.exists(example_pdb):
        pytest.skip("Example protein file not found")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Find cavities
        cavities, cavity_data = find_cavities(
            example_pdb,
            output_dir=tmpdir,
        )
        
        if len(cavities) == 0:
            pytest.skip("No cavities found in test protein")
        
        # Fill with water using Monte Carlo
        output_file = fill_cavities_with_water(
            example_pdb,
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
