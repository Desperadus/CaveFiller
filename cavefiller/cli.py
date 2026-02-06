"""Command-line interface for CaveFiller using Typer."""

import typer
from pathlib import Path
from typing import Optional, List
from cavefiller.cavity_finder import find_cavities
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
    probe_in: float = typer.Option(
        1.4,
        help="Probe In radius for cavity detection (Å)",
    ),
    probe_out: float = typer.Option(
        4.0,
        help="Probe Out radius for cavity detection (Å)",
    ),
    volume_cutoff: float = typer.Option(
        5.0,
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
):
    """
    Find cavities in a protein and fill them with water molecules.
    
    This tool performs the following steps:
    1. Detects cavities in the protein using KVFinder
    2. Allows user to select which cavities to fill
    3. Uses Packmol to fill selected cavities with water molecules
    """
    typer.echo(f"🔍 Analyzing protein: {protein_file}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Find cavities using pykvfinder
    typer.echo("Step 1: Finding cavities with KVFinder...")
    cavities, cavity_data = find_cavities(
        str(protein_file),
        probe_in=probe_in,
        probe_out=probe_out,
        volume_cutoff=volume_cutoff,
        output_dir=str(output_dir),
    )
    
    if not cavities:
        typer.echo("❌ No cavities found in the protein.", err=True)
        raise typer.Exit(code=1)
    
    typer.echo(f"✅ Found {len(cavities)} cavities")
    
    # Step 2: Select cavities to fill
    typer.echo("\nStep 2: Selecting cavities to fill...")
    
    if cavity_ids:
        # Parse cavity IDs from command line
        selected_ids = [int(x.strip()) for x in cavity_ids.split(",")]
        selected_cavities = [c for c in cavities if c["id"] in selected_ids]
        if not selected_cavities:
            typer.echo(f"❌ No cavities found with IDs: {cavity_ids}", err=True)
            raise typer.Exit(code=1)
    elif auto_select:
        # Auto-select all cavities
        selected_cavities = cavities
        typer.echo(f"Auto-selecting all {len(cavities)} cavities")
    else:
        # Interactive selection
        selected_cavities = select_cavities(cavities)
    
    if not selected_cavities:
        typer.echo("❌ No cavities selected.", err=True)
        raise typer.Exit(code=1)
    
    typer.echo(f"✅ Selected {len(selected_cavities)} cavities")
    
    # Step 3: Fill cavities with water
    typer.echo("\nStep 3: Filling cavities with water using Packmol...")
    output_file = fill_cavities_with_water(
        str(protein_file),
        selected_cavities,
        cavity_data,
        str(output_dir),
    )
    
    typer.echo(f"✅ Success! Output saved to: {output_file}")
    typer.echo("\n🎉 CaveFiller completed successfully!")


if __name__ == "__main__":
    app()
