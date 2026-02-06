"""Interactive cavity selection interface."""

from typing import List, Dict, Any
import sys


def select_cavities(cavities: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Allow user to select which cavities to fill.
    
    Args:
        cavities: List of cavity dictionaries with id, volume, and area
        
    Returns:
        List of selected cavity dictionaries
    """
    print("\n" + "=" * 60)
    print("Available Cavities:")
    print("=" * 60)
    print(f"{'ID':<6} {'Volume (Ų)':<15} {'Area (ų)':<15}")
    print("-" * 60)
    
    for cavity in cavities:
        print(
            f"{cavity['id']:<6} {cavity['volume']:<15.2f} {cavity['area']:<15.2f}"
        )
    
    print("=" * 60)
    print("\nEnter cavity IDs to fill (comma-separated, e.g., '1,2,3')")
    print("Or enter 'all' to select all cavities")
    print("Or enter 'q' to quit")
    
    while True:
        try:
            user_input = input("\nYour selection: ").strip().lower()
            
            if user_input == 'q':
                print("Selection cancelled.")
                return []
            
            if user_input == 'all':
                print(f"Selected all {len(cavities)} cavities")
                return cavities
            
            # Parse comma-separated IDs
            selected_ids = [int(x.strip()) for x in user_input.split(",")]
            
            # Validate IDs
            valid_ids = {c["id"] for c in cavities}
            invalid_ids = [sid for sid in selected_ids if sid not in valid_ids]
            
            if invalid_ids:
                print(f"Invalid cavity IDs: {invalid_ids}")
                print(f"Valid IDs are: {sorted(valid_ids)}")
                continue
            
            # Get selected cavities
            selected = [c for c in cavities if c["id"] in selected_ids]
            
            if selected:
                print(f"\nSelected {len(selected)} cavities: {[c['id'] for c in selected]}")
                return selected
            else:
                print("No cavities selected. Please try again.")
                
        except ValueError as e:
            print(f"Invalid input: {e}")
            print("Please enter comma-separated numbers, 'all', or 'q'")
        except EOFError:
            # Handle case when stdin is not available (e.g., in tests)
            print("\nNo input available. Selecting all cavities by default.")
            return cavities
        except KeyboardInterrupt:
            print("\n\nSelection cancelled.")
            return []
