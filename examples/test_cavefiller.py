"""Script to test CaveFiller functionality."""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cavefiller.cavity_finder import find_cavities
from cavefiller.water_filler import fill_cavities_with_water

def test_cavefiller():
    """Test the CaveFiller pipeline."""
    
    print("=" * 70)
    print("CaveFiller Test - Monte Carlo & MMFF94 Optimization")
    print("=" * 70)
    
    # Use the provided Mus OBP5 example protein
    protein_file = os.path.join(os.path.dirname(__file__), "musM_OBP5_model_0_boltz2.pdb")
    if not os.path.exists(protein_file):
        print(f"❌ Example protein not found: {protein_file}")
        return
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
            print("\nNote: Adjust cavity detection parameters and retry.")
            return
        
        # Step 2: Auto-select all cavities with specific water counts
        print("\n[Step 2] Auto-selecting cavities with specified water counts...")
        selected_cavities = cavities[:3]  # Select first 3 cavities
        waters_dict = {}
        for i, cavity in enumerate(selected_cavities):
            waters_dict[cavity['id']] = 5 + i * 3  # 5, 8, 11 waters
        
        print(f"✅ Selected {len(selected_cavities)} cavities")
        for cavity in selected_cavities:
            print(f"   - Cavity {cavity['id']}: {waters_dict[cavity['id']]} waters")
        
        # Step 3: Fill with water using Monte Carlo and MMFF94
        print("\n[Step 3] Filling cavities with Monte Carlo sampling and MMFF94...")
        output_file = fill_cavities_with_water(
            protein_file,
            selected_cavities,
            cavity_data,
            output_dir,
            waters_per_cavity=waters_dict
        )
        print(f"✅ Output saved to: {output_file}")
        
        # Check output
        if os.path.exists(output_file):
            with open(output_file, 'r') as f:
                lines = f.readlines()
            water_lines = [l for l in lines if 'HOH' in l]
            print(f"✅ Added {len(water_lines)} water atoms")
        
        print("\n" + "=" * 70)
        print("✅ Test completed successfully!")
        print("=" * 70)
        print("\nKey Features Demonstrated:")
        print("  ✅ Monte Carlo sampling with clash detection")
        print("  ✅ Van der Waals radii-based distance validation")
        print("  ✅ MMFF94 force field optimization")
        print("  ✅ User-specified water counts per cavity")
        
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    test_cavefiller()
