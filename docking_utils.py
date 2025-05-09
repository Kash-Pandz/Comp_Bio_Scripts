import subprocess
import numpy as np
import MDAnalysis as mda


def split_pdb(pdb_file, ligand_name, protein_file="protein.pdb", ligand_file="ligand.pdb"):
    """Splits input PDB file into protein and ligand files."""

    # Load PDB file into MDAnalysis universe
    universe = mda.Universe(pdb_file)

    # Select protein and ligand atoms
    protein = universe.select_atoms("protein")
    ligand = universe.select_atoms(f"resname {ligand_name}")

    # Write protein and ligand pdb files
    protein.write(protein_file)
    ligand.write(ligand_file)


def protein_to_pdbqt(prot_pdb, prot_pdbqt):
   """Convert protein pdb into a pdbqt file."""
    command = [
        'obabel', 
        '-ipdb', prot_pdb,
        '-h',
        '--partialcharge', 'eem',
        '-opdbqt', 
        '-O', prot_pdbqt, 
        '-xr'
    ]
    subprocess.run(command, stdout=subprocess.DEVNULL)


def ligand_to_pdbqt(ligand_pdb, ligand_pdbqt):
   """Convert ligand sdf or pdb file into pdbqt format."""
    command = [
        'obabel', 
        '-ipdb', lig_pdb,
        '-h',
        '--partialcharge', 'eem',
        '-opdbqt', 
        '-O', lig_pdbqt, 
    ]
    subprocess.run(command, stdout=subprocess.DEVNULL)

