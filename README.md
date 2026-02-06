# CaveFiller

A Python tool to find and fill protein cavities with water molecules using KVFinder and Packmol.

## Features

- 🔍 **Cavity Detection**: Uses pyKVFinder to detect cavities in protein structures
- 🎯 **Interactive Selection**: Select specific cavities to fill or auto-select all
- 💧 **Water Filling**: Fills selected cavities with water molecules using Packmol
- 🖥️ **CLI Interface**: Easy-to-use command-line interface built with Typer

## Installation

### Prerequisites

1. **Python**: Python 3.8 or higher
2. **Packmol** (optional but recommended): Install from [Packmol website](http://m3g.iqm.unicamp.br/packmol/home.shtml)

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
4. Fill the selected cavities with water molecules
5. Save the output to `./output/protein_filled.pdb`

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
- `--volume-cutoff FLOAT`: Minimum cavity volume to consider in ų (default: 5.0)
- `--auto-select`: Automatically select all cavities without user interaction
- `--cavity-ids TEXT`: Comma-separated list of cavity IDs to fill (e.g., '1,2,3')

### Examples

**Interactive cavity selection:**
```bash
cavefiller protein.pdb --output-dir results
```

**Auto-select all cavities:**
```bash
cavefiller protein.pdb --auto-select
```

**Fill specific cavities:**
```bash
cavefiller protein.pdb --cavity-ids "1,3,5"
```

**Custom cavity detection parameters:**
```bash
cavefiller protein.pdb --probe-in 1.2 --probe-out 5.0 --volume-cutoff 10.0
```

## Workflow

1. **Cavity Detection**: The tool uses pyKVFinder to detect cavities in the input protein structure
2. **Cavity Analysis**: Displays information about detected cavities (ID, volume, surface area)
3. **Cavity Selection**: 
   - Interactive mode: User selects cavities by entering IDs
   - Auto mode: All cavities are selected automatically
   - Command-line mode: Specific cavities are pre-selected
4. **Water Filling**: 
   - If Packmol is available: Uses Packmol to optimally pack water molecules
   - Fallback mode: Uses simple grid-based water placement

## Output

The tool generates the following files in the output directory:

- `protein_filled.pdb`: Protein structure with water molecules in selected cavities
- `cavities.toml`: Detailed cavity detection results from KVFinder
- `packmol.inp`: Packmol input file (if using Packmol)
- `water.pdb`: Water molecule template (if using Packmol)

## Dependencies

- **typer**: CLI framework
- **pyKVFinder**: Cavity detection
- **rdkit**: Molecular manipulation
- **numpy**: Numerical operations
- **biopython**: PDB file handling
- **Packmol** (optional): Optimal water molecule packing

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
- Packmol: Martínez et al. (2009) J. Comput. Chem.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.