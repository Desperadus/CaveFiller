"""Command-line interface for CaveFiller using Typer."""

import typer
from pathlib import Path
from typing import Optional
from cavefiller.cavity_finder import (
    find_cavities,
    DEFAULT_GRID_STEP,
    DEFAULT_PROBE_IN,
    DEFAULT_PROBE_OUT,
    DEFAULT_EXTERIOR_TRIM_DISTANCE,
    DEFAULT_VOLUME_CUTOFF,
)
from cavefiller.cavity_selector import select_cavities
from cavefiller.water_filler import fill_cavities_with_water

app = typer.Typer(help="CaveFiller - Find and fill protein cavities with water molecules")


@app.command()
def run(
    protein_file: Path = typer.Argument(
        ...,
        exists=True,
        help="Path to the protein PDB file",
    ),
    output_dir: Path = typer.Option(
        Path("./output"),
        help="Directory to save output files",
    ),
    grid_step: float = typer.Option(
        DEFAULT_GRID_STEP,
        "--grid-step",
        help="Grid spacing for cavity detection (Å)",
    ),
    probe_in: float = typer.Option(
        DEFAULT_PROBE_IN,
        help="Probe In radius for cavity detection (Å)",
    ),
    probe_out: float = typer.Option(
        DEFAULT_PROBE_OUT,
        help="Probe Out radius for cavity detection (Å)",
    ),
    exterior_trim_distance: float = typer.Option(
        DEFAULT_EXTERIOR_TRIM_DISTANCE,
        "--exterior-trim-distance",
        help="Exterior trim distance (KVFinder removal distance) (Å)",
    ),
    volume_cutoff: float = typer.Option(
        DEFAULT_VOLUME_CUTOFF,
        help="Minimum cavity volume to consider (Å³)",
    ),
    auto_select: bool = typer.Option(
        False,
        help="Automatically select all cavities (no user interaction)",
    ),
    cavity_ids: Optional[str] = typer.Option(
        None,
        help="Comma-separated list of cavity IDs to fill (e.g., '1,2,3'). If not provided, user will be prompted.",
    ),
    waters_per_cavity: Optional[str] = typer.Option(
        None,
        help="Comma-separated list of water counts per cavity (e.g., '10,15,20'). Must match cavity_ids order.",
    ),
    optimize_mmff94: bool = typer.Option(
        True,
        "--optimize-mmff94/--no-optimize-mmff94",
        help="Run MMFF94 after placement with protein atoms fixed and waters movable.",
    ),
    mmff_max_iterations: int = typer.Option(
        300,
        help="Maximum MMFF94 iterations when optimization is enabled.",
    ),
):
    """
    Find cavities in a protein and fill them with explicit water molecules.

    This tool performs the following steps:
    1. Detects cavities in the protein using KVFinder
    2. Allows user to select which cavities to fill
    3. Uses cavity-grid Monte Carlo sampling to place water oxygens
    4. Optionally runs MMFF94 with protein fixed and waters movable
    5. Builds explicit RDKit H-O-H waters and writes a combined PDB
    """
    typer.echo(f"🔍 Analyzing protein: {protein_file}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Find cavities using pykvfinder
    typer.echo("Step 1: Finding cavities with KVFinder...")
    cavities, cavity_data = find_cavities(
        str(protein_file),
        step=grid_step,
        probe_in=probe_in,
        probe_out=probe_out,
        removal_distance=exterior_trim_distance,
        volume_cutoff=volume_cutoff,
        output_dir=str(output_dir),
    )
    
    if not cavities:
        typer.echo("❌ No cavities found in the protein.", err=True)
        raise typer.Exit(code=1)
    
    typer.echo(f"✅ Found {len(cavities)} cavities")
    
    # Step 2: Select cavities to fill
    typer.echo("\nStep 2: Selecting cavities to fill...")
    
    waters_dict = {}
    
    if cavity_ids:
        # Parse cavity IDs from command line
        selected_ids = [int(x.strip()) for x in cavity_ids.split(",")]
        selected_cavities = [c for c in cavities if c["id"] in selected_ids]
        if not selected_cavities:
            typer.echo(f"❌ No cavities found with IDs: {cavity_ids}", err=True)
            raise typer.Exit(code=1)
        
        # Parse waters per cavity if provided
        if waters_per_cavity:
            water_counts = [int(x.strip()) for x in waters_per_cavity.split(",")]
            if len(water_counts) != len(selected_ids):
                typer.echo("❌ Number of water counts must match number of cavity IDs", err=True)
                raise typer.Exit(code=1)
            waters_dict = dict(zip(selected_ids, water_counts))
            
    elif auto_select:
        # Auto-select all cavities
        selected_cavities = cavities
        typer.echo(f"Auto-selecting all {len(cavities)} cavities")
        # Use default water counts
        for cavity in selected_cavities:
            waters_dict[cavity['id']] = max(1, int(cavity['volume'] / 30))
    else:
        # Interactive selection
        selected_cavities, waters_dict = select_cavities(cavities, prompt_for_waters=True)
    
    if not selected_cavities:
        typer.echo("❌ No cavities selected.", err=True)
        raise typer.Exit(code=1)
    
    typer.echo(f"✅ Selected {len(selected_cavities)} cavities")
    
    # Step 3: Fill cavities with water
    typer.echo("\nStep 3: Filling cavities with water using Monte Carlo sampling...")
    typer.echo("         (with clash detection, optional fixed-protein MMFF94, and RDKit HOH generation)")
    output_file = fill_cavities_with_water(
        str(protein_file),
        selected_cavities,
        cavity_data,
        str(output_dir),
        waters_per_cavity=waters_dict,
        optimize_mmff94=optimize_mmff94,
        mmff_max_iterations=mmff_max_iterations,
    )
    
    typer.echo(f"\n✅ Success! Output saved to: {output_file}")
    typer.echo("\n🎉 CaveFiller completed successfully!")


if __name__ == "__main__":
    app()
