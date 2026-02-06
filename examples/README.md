# CaveFiller Examples

This directory contains example files and scripts to demonstrate CaveFiller functionality.

## Files

- `musM_OBP5_model_0_boltz2.pdb`: Example Mus OBP5 structure for cavity filling tests
- `test_cavefiller.py`: A comprehensive test script that demonstrates the full workflow

## Running the Test Script

```bash
python examples/test_cavefiller.py
```

This will:
1. Load the provided Mus OBP5 example protein
2. Run cavity detection
3. Select cavities and set water counts
4. Fill cavities with water molecules
5. Display results

## Using the CLI

### Example 1: Interactive Mode

```bash
cavefiller examples/musM_OBP5_model_0_boltz2.pdb --output-dir examples/output
```

This will prompt you to select which cavities to fill.

### Example 2: Auto-select All Cavities (not generally recommended)

```bash
cavefiller examples/musM_OBP5_model_0_boltz2.pdb --auto-select --output-dir examples/output
```

### Example 3: Select Specific Cavities

```bash
cavefiller examples/musM_OBP5_model_0_boltz2.pdb --cavity-ids "1,2,3" --output-dir examples/output
```

### Example 4: Custom Detection Parameters

```bash
cavefiller examples/musM_OBP5_model_0_boltz2.pdb \
  --probe-in 1.2 \
  --probe-out 5.0 \
  --volume-cutoff 10.0 \
  --cavity-ids "1,2"
```

### Example 5: Recommended Manual Selection + MMFF94

```bash
cavefiller examples/musM_OBP5_model_0_boltz2.pdb --output-dir examples/output
```

Notes:
- Prefer manual cavity selection and manual water counts. `--auto-select` often places too many waters.
- Keep MMFF94 optimization enabled (`--optimize-mmff94`, default) for better final geometries.

## Expected Output

The tool will create the following files in the output directory:

- `protein_filled.pdb`: The protein structure with water molecules added to selected cavities
- Other intermediate files depending on the method used (Packmol or simple placement)

## Notes

- Real protein structures from PDB can be downloaded and used with this tool
