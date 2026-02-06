# CaveFiller

A Python tool to find and fill protein cavities with water molecules using KVFinder, Monte Carlo sampling, and MMFF94 optimization.

## Features

- 🔍 **Cavity Detection**: Uses pyKVFinder to detect cavities in protein structures
- 🎯 **Interactive Selection**: Select specific cavities to fill or auto-select all
- 🎲 **Monte Carlo Sampling**: Places water molecules using Monte Carlo sampling with clash detection
- ⚛️ **MMFF94 Optimization**: Optimizes water positions using RDKit's MMFF94 force field
- 🖥️ **CLI Interface**: Easy-to-use command-line interface built with Typer

## Installation

### Prerequisites

1. **Python**: Python 3.8 or higher

### Install CaveFiller

```bash
# Clone the repository
git clone https://github.com/Desperadus/CaveFiller.git
cd CaveFiller

# Install the package
pip install -e .
```

## Usage

### Basic Usage

```bash
cavefiller protein.pdb
```

This will:
1. Detect cavities in `protein.pdb`
2. Display a list of found cavities with their volumes and areas
3. Prompt you to select which cavities to fill
4. Prompt you for the number of water molecules per cavity
5. Place waters using Monte Carlo sampling with clash detection
6. Optimize water positions using MMFF94 force field
7. Save the output to `./output/protein_filled.pdb`

### Command-line Options

```bash
cavefiller [PROTEIN_FILE] [OPTIONS]
```

**Arguments:**
- `PROTEIN_FILE`: Path to the protein PDB file (required)

**Options:**
- `--output-dir PATH`: Directory to save output files (default: `./output`)
- `--probe-in FLOAT`: Probe In radius for cavity detection in Ångströms (default: 1.4)
- `--probe-out FLOAT`: Probe Out radius for cavity detection in Ångströms (default: 4.0)
- `--volume-cutoff FLOAT`: Minimum cavity volume to consider in Ų (default: 5.0)
- `--auto-select`: Automatically select all cavities without user interaction
- `--cavity-ids TEXT`: Comma-separated list of cavity IDs to fill (e.g., '1,2,3')
- `--waters-per-cavity TEXT`: Comma-separated list of water counts (e.g., '10,15,20'), must match cavity-ids order

### Examples

**Interactive cavity and water selection:**
```bash
cavefiller protein.pdb --output-dir results
```

**Auto-select all cavities with default water counts:**
```bash
cavefiller protein.pdb --auto-select
```

**Fill specific cavities with specific water counts:**
```bash
cavefiller protein.pdb --cavity-ids "1,3,5" --waters-per-cavity "10,15,20"
```

**Custom cavity detection parameters:**
```bash
cavefiller protein.pdb --probe-in 1.2 --probe-out 5.0 --volume-cutoff 10.0
```

## Workflow

1. **Cavity Detection**: The tool uses pyKVFinder to detect cavities in the input protein structure
2. **Cavity Analysis**: Displays information about detected cavities (ID, volume, surface area)
3. **Cavity Selection**: 
   - Interactive mode: User selects cavities and specifies water counts
   - Auto mode: All cavities are selected with automatic water count estimation
   - Command-line mode: Specific cavities and water counts are pre-selected
4. **Water Placement**: 
   - Monte Carlo sampling places waters randomly in cavity
   - Clash detection validates each position against protein atoms and other waters
   - Uses Van der Waals radii for distance calculations
5. **MMFF94 Optimization**:
   - Waters are optimized using RDKit's MMFF94 force field
   - Ensures reasonable geometry and hydrogen bonding
   - Waters that drift outside cavity are filtered out

## Algorithm Details

### Monte Carlo Sampling
- Randomly samples positions within cavity bounding box
- Validates position is within cavity (< 1.5 Å from cavity grid point)
- Checks for clashes with protein atoms (minimum distance based on VDW radii)
- Checks for clashes with other waters (minimum 2.8 Å separation)
- Attempts up to 1000 placements per water molecule

### Clash Detection
- Uses Van der Waals radii for different atom types (H, C, N, O, S, P)
- Minimum water-protein distance: 2.4 Å
- Minimum water-water distance: 2.8 Å
- Tolerance of 0.5 Å for VDW overlap

### MMFF94 Optimization
- Creates proper H-O-H geometry for each water
- Optimizes all waters simultaneously
- Maximum 200 iterations per optimization
- Filters waters that move > 2.0 Å from cavity

## Output

The tool generates the following files in the output directory:

- `protein_filled.pdb`: Protein structure with optimized water molecules in selected cavities

## Dependencies

- **typer**: CLI framework
- **pyKVFinder**: Cavity detection
- **rdkit**: Molecular manipulation and MMFF94 optimization
- **numpy**: Numerical operations
- **biopython**: PDB file handling

## Development

### Running Tests

```bash
pip install -e ".[dev]"
pytest
```

### Code Formatting

```bash
black cavefiller/
ruff check cavefiller/
```

## License

See LICENSE file for details.

## Citation

If you use CaveFiller in your research, please cite:

- pyKVFinder: Guerra et al. (2020) BMC Bioinformatics
- RDKit: RDKit: Open-source cheminformatics; http://www.rdkit.org

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.