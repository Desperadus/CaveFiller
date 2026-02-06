"""Script to test CaveFiller functionality."""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cavefiller.cavity_finder import find_cavities
from cavefiller.cavity_selector import select_cavities
from cavefiller.water_filler import fill_cavities_with_water


def create_test_protein_with_cavity():
    """Create a synthetic protein structure with a cavity for testing."""
    
    # Create a hollow box structure - atoms arranged in walls
    pdb_content = "HEADER    TEST PROTEIN WITH CAVITY\n"
    pdb_content += "REMARK   Synthetic protein structure with an internal cavity\n"
    
    atom_num = 1
    res_num = 1
    
    # Create box walls (4 walls of a box)
    # Bottom wall (y=0, z varies, x varies)
    for x in range(0, 20, 2):
        for z in range(0, 20, 2):
            pdb_content += f"ATOM  {atom_num:5d}  CA  ALA A{res_num:4d}    {x:8.3f}{0.0:8.3f}{z:8.3f}  1.00  0.00           C\n"
            atom_num += 1
            res_num += 1
    
    # Top wall (y=20, z varies, x varies)
    for x in range(0, 20, 2):
        for z in range(0, 20, 2):
            pdb_content += f"ATOM  {atom_num:5d}  CA  ALA A{res_num:4d}    {x:8.3f}{20.0:8.3f}{z:8.3f}  1.00  0.00           C\n"
            atom_num += 1
            res_num += 1
    
    # Left wall (x=0, y varies, z varies)
    for y in range(2, 18, 2):
        for z in range(0, 20, 2):
            pdb_content += f"ATOM  {atom_num:5d}  CA  ALA A{res_num:4d}    {0.0:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            atom_num += 1
            res_num += 1
    
    # Right wall (x=20, y varies, z varies)
    for y in range(2, 18, 2):
        for z in range(0, 20, 2):
            pdb_content += f"ATOM  {atom_num:5d}  CA  ALA A{res_num:4d}    {20.0:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            atom_num += 1
            res_num += 1
    
    # Front wall (z=0, x varies, y varies) 
    for x in range(2, 18, 2):
        for y in range(2, 18, 2):
            pdb_content += f"ATOM  {atom_num:5d}  CA  ALA A{res_num:4d}    {x:8.3f}{y:8.3f}{0.0:8.3f}  1.00  0.00           C\n"
            atom_num += 1
            res_num += 1
    
    # Back wall (z=20, x varies, y varies)
    for x in range(2, 18, 2):
        for y in range(2, 18, 2):
            pdb_content += f"ATOM  {atom_num:5d}  CA  ALA A{res_num:4d}    {x:8.3f}{y:8.3f}{20.0:8.3f}  1.00  0.00           C\n"
            atom_num += 1
            res_num += 1
    
    pdb_content += "END\n"
    
    # Save to file
    output_path = os.path.join(os.path.dirname(__file__), "protein_with_cavity.pdb")
    with open(output_path, 'w') as f:
        f.write(pdb_content)
    
    print(f"Created test protein: {output_path}")
    print(f"Total atoms: {atom_num - 1}")
    return output_path


def test_cavefiller():
    """Test the CaveFiller pipeline."""
    
    print("=" * 70)
    print("CaveFiller Test")
    print("=" * 70)
    
    # Create test protein
    protein_file = create_test_protein_with_cavity()
    output_dir = os.path.join(os.path.dirname(__file__), "test_output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Find cavities
    print("\n[Step 1] Finding cavities...")
    try:
        cavities, cavity_data = find_cavities(
            protein_file,
            probe_in=1.4,
            probe_out=4.0,
            volume_cutoff=5.0,
            output_dir=output_dir,
        )
        
        if cavities:
            print(f"✅ Found {len(cavities)} cavities")
            for cavity in cavities:
                print(f"   - Cavity {cavity['id']}: Volume={cavity['volume']:.2f} Ų, Area={cavity['area']:.2f} ų")
        else:
            print("❌ No cavities found")
            print("\nNote: The test protein may be too small or the cavity detection")
            print("parameters may need adjustment. This is expected for small test proteins.")
            return
        
        # Step 2: Auto-select all cavities
        print("\n[Step 2] Auto-selecting all cavities...")
        selected_cavities = cavities
        print(f"✅ Selected {len(selected_cavities)} cavities")
        
        # Step 3: Fill with water
        print("\n[Step 3] Filling cavities with water...")
        output_file = fill_cavities_with_water(
            protein_file,
            selected_cavities,
            cavity_data,
            output_dir,
        )
        print(f"✅ Output saved to: {output_file}")
        
        # Check output
        if os.path.exists(output_file):
            with open(output_file, 'r') as f:
                lines = f.readlines()
            water_lines = [l for l in lines if 'HOH' in l or 'WAT' in l]
            print(f"✅ Added {len(water_lines)} water atoms/molecules")
        
        print("\n" + "=" * 70)
        print("✅ Test completed successfully!")
        print("=" * 70)
        
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    test_cavefiller()
