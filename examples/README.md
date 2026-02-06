# CaveFiller Examples

This directory contains example files and scripts to demonstrate CaveFiller functionality.

## Files

- `protein_with_cavity.pdb`: A synthetic protein structure with an internal cavity, ideal for testing
- `sample_protein.pdb`: A small sample protein structure (may be too small for cavity detection)
- `test_cavefiller.py`: A comprehensive test script that demonstrates the full workflow

## Running the Test Script

```bash
python examples/test_cavefiller.py
```

This will:
1. Create a test protein with a cavity
2. Run cavity detection
3. Auto-select all cavities
4. Fill cavities with water molecules
5. Display results

## Using the CLI

### Example 1: Interactive Mode

```bash
cavefiller examples/protein_with_cavity.pdb --output-dir examples/output
```

This will prompt you to select which cavities to fill.

### Example 2: Auto-select All Cavities

```bash
cavefiller examples/protein_with_cavity.pdb --auto-select --output-dir examples/output
```

### Example 3: Select Specific Cavities

```bash
cavefiller examples/protein_with_cavity.pdb --cavity-ids "1,2,3" --output-dir examples/output
```

### Example 4: Custom Detection Parameters

```bash
cavefiller examples/protein_with_cavity.pdb \
  --probe-in 1.2 \
  --probe-out 5.0 \
  --volume-cutoff 10.0 \
  --auto-select
```

## Expected Output

The tool will create the following files in the output directory:

- `protein_filled.pdb`: The protein structure with water molecules added to selected cavities
- Other intermediate files depending on the method used (Packmol or simple placement)

## Notes

- If Packmol is not installed, the tool will fall back to a simple grid-based water placement method
- The synthetic test protein creates a hollow box structure that reliably has detectable cavities
- Real protein structures from PDB can be downloaded and used with this tool
