"""Unit tests for CaveFiller."""

import pytest
import os
import tempfile
from pathlib import Path


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


def test_water_template_creation():
    """Test that water template can be created."""
    from cavefiller.water_filler import create_water_template
    
    with tempfile.TemporaryDirectory() as tmpdir:
        water_file = create_water_template(tmpdir)
        assert os.path.exists(water_file)
        
        with open(water_file, 'r') as f:
            content = f.read()
            assert 'HOH' in content
            assert 'HETATM' in content


def test_packmol_availability_check():
    """Test packmol availability check."""
    from cavefiller.water_filler import is_packmol_available
    
    # This will be False in most test environments
    result = is_packmol_available()
    assert isinstance(result, bool)


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


def test_simple_water_filling():
    """Test simple water filling method."""
    from cavefiller.water_filler import fill_with_simple_placement
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
        
        # Fill with water
        output_file = os.path.join(tmpdir, "filled.pdb")
        result = fill_with_simple_placement(
            example_pdb,
            cavities[:1],  # Just first cavity
            cavity_data,
            output_file,
        )
        
        assert os.path.exists(result)
        
        # Check that water was added
        with open(result, 'r') as f:
            content = f.read()
            assert 'HOH' in content or 'WAT' in content
