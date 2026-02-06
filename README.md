# CaveFiller

A Python tool to find and fill protein cavities with water molecules using KVFinder, Monte Carlo sampling, and RDKit-based explicit water generation.

## Features

-  **Cavity Detection**: Uses pyKVFinder to detect cavities in protein structures
-  **Interactive Selection**: Select specific cavities to fill with user-defined water counts
-  **Monte Carlo Sampling**: Places water molecules using Monte Carlo sampling with clash detection
-  **Explicit Waters**: Builds full H-O-H waters with RDKit (including hydrogens)
-  **CLI Interface**: Easy-to-use command-line interface built with Typer

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
6. Build explicit RDKit H-O-H waters and export a combined PDB
7. Save the output to `./output/protein_filled.pdb`

### Command-line Options

```bash
cavefiller [PROTEIN_FILE] [OPTIONS]
```

**Arguments:**
- `PROTEIN_FILE`: Path to the protein PDB file (required)

**Options:**
- `--output-dir PATH`: Directory to save output files (default: `./output`)
- `--grid-step FLOAT`: Grid spacing for cavity detection in Ångströms (default: 0.6)
- `--probe-in FLOAT`: Probe In radius for cavity detection in Ångströms (default: 1.4)
- `--probe-out FLOAT`: Probe Out radius for cavity detection in Ångströms (default: 4.0)
- `--exterior-trim-distance FLOAT`: Exterior trim distance in Ångströms (default: 2.4)
- `--volume-cutoff FLOAT`: Minimum cavity volume to consider in Ų (default: 5.0)
- `--auto-select`: Automatically select all cavities without user interaction
- `--cavity-ids TEXT`: Comma-separated list of cavity IDs to fill (e.g., '1,2,3')
- `--waters-per-cavity TEXT`: Comma-separated list of water counts (e.g., '10,15,20'), must match cavity-ids order
- `--optimize-mmff94 / --no-optimize-mmff94`: Enable/disable MMFF94 with protein fixed (default: enabled)
- `--mmff-max-iterations INTEGER`: Max MMFF94 iterations (default: 300)
- `--remove-after-optim / --no-remove-after-optim`: After MMFF94, remove waters that fail post-checks (default: enabled)
  - Also accepted: `--remove_after_optim / --no_remove_after_optim`

Recommended usage:
- Prefer interactive/manual cavity and water-count selection over `--auto-select`. Auto-selection often overfills cavities with too many waters.
- Keep `--optimize-mmff94` enabled (recommended) to refine water placement after Monte Carlo sampling.
- Use `--no-remove-after-optim` if you want to keep all waters after MMFF94, even if they clash or move out of cavity bounds.

### Examples

**Interactive cavity and water selection:**
```bash
cavefiller protein.pdb --output-dir results
```

**Auto-select all cavities with default water counts (not generally recommended):**
```bash
cavefiller protein.pdb --auto-select
```

**Fill specific cavities with specific water counts:**
```bash
cavefiller protein.pdb --cavity-ids "1,3,5" --waters-per-cavity "10,15,20"
```

**Custom cavity detection parameters:**
```bash
cavefiller protein.pdb --grid-step 0.6 --probe-in 1.4 --probe-out 4.0 --exterior-trim-distance 2.4 --volume-cutoff 5.0
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
5. **RDKit Water Construction**:
   - Explicit H-O-H waters are generated with ideal geometry
   - Waters include hydrogens and proper HOH residue records in the output PDB

## Algorithm Details

### Monte Carlo Sampling
- Samples around cavity grid points with small local jitter
- Validates position stays near cavity voxels (< 0.7 Å from a grid point)
- Checks for clashes with protein atoms (minimum distance based on VDW radii)
- Checks for clashes with other waters (minimum 2.7 Å separation)
- Attempts up to 500 placements per water molecule

### Clash Detection
- Uses Van der Waals radii for different atom types (H, C, N, O, S, P)
- Minimum water-protein distance: 2.35 Å
- Minimum water-water distance: 2.7 Å
- Tolerance of 0.5 Å for VDW overlap

### RDKit Water Geometry
- Creates proper H-O-H geometry for each water
- Writes explicit HOH residues (O, H1, H2) into output PDB

## Output

The tool generates the following files in the output directory:

- `protein_filled.pdb`: Protein structure with explicit water molecules in selected cavities

## Dependencies

- **typer**: CLI framework
- **pyKVFinder**: Cavity detection
- **rdkit**: Molecular manipulation and explicit water generation
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

### Automated CI/CD and PyPI Publishing

This repository includes GitHub Actions workflow at `.github/workflows/ci-cd.yml` that:
- Runs `pytest` on every push to `main`
- Runs `pytest` on every pull request targeting `main`
- Builds package distributions after tests pass
- Publishes to PyPI only on pushes to `main` where `pyproject.toml` `project.version` changed

#### One-time setup for automatic PyPI publishing

1. Create a PyPI account at https://pypi.org and create your project once (or publish once manually so the name exists).
2. In PyPI, open your project settings and add a **Trusted Publisher**:
   - Owner: your GitHub username/org
   - Repository: `Desperadus/CaveFiller`
   - Workflow name: `CI/CD`
   - Environment: leave empty (unless you choose to use one)
3. In GitHub, ensure Actions are enabled for the repository.

No PyPI API token secret is needed when using Trusted Publishing.

#### Releasing a new version

1. Bump version in both:
   - `pyproject.toml` (`project.version`)
   - `cavefiller/__init__.py` (`__version__`)
2. Commit and push to `main`.
3. CI will publish that pushed version to PyPI automatically, but only if `pyproject.toml` version changed versus the previous commit on `main`.

## License

See LICENSE file for details.

## Citation

If you use CaveFiller in your research, please cite:

- pyKVFinder: Guerra et al. (2020) BMC Bioinformatics
- RDKit: RDKit: Open-source cheminformatics; http://www.rdkit.org

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
