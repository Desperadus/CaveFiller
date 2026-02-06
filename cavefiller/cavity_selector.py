"""Interactive cavity selection interface."""

from typing import List, Dict, Any, Optional
import sys


def select_cavities(cavities: List[Dict[str, Any]], prompt_for_waters: bool = True) -> tuple:
    """
    Allow user to select which cavities to fill and how many waters per cavity.
    
    Args:
        cavities: List of cavity dictionaries with id, volume, and area
        prompt_for_waters: Whether to prompt for number of waters per cavity
        
    Returns:
        Tuple of (selected cavity dictionaries, waters_per_cavity dict)
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
    
    selected_cavities = None
    while selected_cavities is None:
        try:
            user_input = input("\nYour selection: ").strip().lower()
            
            if user_input == 'q':
                print("Selection cancelled.")
                return [], {}
            
            if user_input == 'all':
                print(f"Selected all {len(cavities)} cavities")
                selected_cavities = cavities
            else:
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
                selected_cavities = [c for c in cavities if c["id"] in selected_ids]
                
                if selected_cavities:
                    print(f"\nSelected {len(selected_cavities)} cavities: {[c['id'] for c in selected_cavities]}")
                else:
                    print("No cavities selected. Please try again.")
                    selected_cavities = None
                    
        except ValueError as e:
            print(f"Invalid input: {e}")
            print("Please enter comma-separated numbers, 'all', or 'q'")
        except EOFError:
            # Handle case when stdin is not available (e.g., in tests)
            print("\nNo input available. Selecting all cavities by default.")
            selected_cavities = cavities
        except KeyboardInterrupt:
            print("\n\nSelection cancelled.")
            return [], {}
    
    # Prompt for number of waters per cavity
    waters_per_cavity = {}
    if prompt_for_waters and selected_cavities:
        print("\n" + "=" * 60)
        print("Specify number of water molecules per cavity")
        print("=" * 60)
        
        for cavity in selected_cavities:
            # Default estimate based on volume
            default_waters = max(1, int(cavity['volume'] / 30))
            
            while True:
                try:
                    prompt = f"Cavity {cavity['id']} (volume: {cavity['volume']:.2f} Ų) - waters [default: {default_waters}]: "
                    user_input = input(prompt).strip()
                    
                    if user_input == '':
                        waters_per_cavity[cavity['id']] = default_waters
                        break
                    else:
                        n_waters = int(user_input)
                        if n_waters < 0:
                            print("  Error: Number of waters must be non-negative")
                            continue
                        waters_per_cavity[cavity['id']] = n_waters
                        break
                        
                except ValueError:
                    print("  Error: Please enter a valid number")
                except EOFError:
                    # Use default
                    waters_per_cavity[cavity['id']] = default_waters
                    break
                except KeyboardInterrupt:
                    print("\n\nCancelled. Using defaults for remaining cavities.")
                    for remaining_cavity in selected_cavities:
                        if remaining_cavity['id'] not in waters_per_cavity:
                            default = max(1, int(remaining_cavity['volume'] / 30))
                            waters_per_cavity[remaining_cavity['id']] = default
                    return selected_cavities, waters_per_cavity
    
    return selected_cavities, waters_per_cavity
