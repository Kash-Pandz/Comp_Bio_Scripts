from pathlib import Path
import os, subprocess, tempfile
import numpy as np
import pandas as pd
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from rdkit import Chem
from rdkit.Chem import AllChem
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.lib.distances import capped_distance
from loguru import logger


def configure_logging(verbose: bool):
    """Configure loguru logging."""
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(lambda msg: print(msg, end=""), level=level)


def run(command, description=""):
    """Run the shell command and stream stdout/stderr in real time."""
    logger.info(f"[RUN] {description}: {command}\n")
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        texts=True
    )
    for line in iter(process.stdout.readline, ''):
        if line:
            logger.info(line.strip())
    process.stdout.close()
    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)
    return return_code


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
    


def extract_center(pdb_file, selection, output_file="centre.txt"):
    """Extracts the geometric center of binding pocket or native ligand from a PDB file."""
    
    # Load PDB file into MDAnalyisis universe
    universe = mda.Universe(pdb_file)

    # Select ligand atoms 
    ligand = universe.select_atoms(f"resname {ligand_name}")

    if len(ligand) > 0:
        cog = ligand.center_of_geometry()
        name = f"Ligand ({selection})"
    else:
        pocket = universe.select_atoms("all")  # Selection made after manual inspection of pocket PDB file obtained from fpocket
        if len(pocket) > 0:
            cog = pocket.center_of_geometry()
            name = "pocket"
        else:
            print("No ligand or pocket found in the PDB file.")
            return None

     return cog.tolist()


### Ligand Preparation ###

def prepare_ligand(smiles_file, output_dir):
    """Convert ligand SMILES to PDBQT files."""
    output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    ligand_pdbqt_files = []

    with open(smiles_file) as f:
        for i, line in enumerate(f):
            fields = line.strip().split()
            if not fields:
                continue

            smi = fields[0]
            lig_name = fields[1] if len(fields) > 1 else f"ligand_{i}"
            cleaned_name = "".join(c if c.isalnum() or c in "_-" else "_" for c in lig_name)
            short_name = (cleaned_name[:3].upper()).ljust(3, "-")

            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                logger.warning(f"Skipping invalid SMILES: {smi}")
                continue
            
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)

            with tempfile.NamedTemporaryFile(suffix=".mol", delete=False) as tmp:
                Chem.MolToMolFile(mol, tmp.name)
                tmp_path = tmp.name

            pdbqt_path = output_dir / f"{short_name}.pdbqt"

            # Convert using Open Babel
            try:
                run(
                    f"obabel {tmp_path} -O {pdbqt_path} --gen3d --partialcharge gasteiger",
                    description=f"Ligand {lig_name}"
                )
            except subprocess.CalledProcessError as e:
                    logger.error(f"Open Babel failed for {lig_name}: {e}")
                    os.remove(tmp_path)
                    continue
            finally:
               # Clean up temporary file
               if os.path.exists(tmp_path):
                    os.remove(tmp_path)

            # Patch residue naem and add REMARK
            with open(pdbqt_path, "r+") as f_out:
                lines = f_out.readlines()
                new_lines = [f"REMARK Name: {short_name}\n"]
                for line in lines:
                    if line.startswith(("HETATM", "ATOM")):
                        line = line[:17] + short_name + line[20:]
                    new_lines.append(line)
                f_out.seek(0)
                f_out.writelines(new_lines)
                f_out.truncate()
            
            ligand_pdbqt_files.append(str(pdbqt_path))
                    
        
    logger.info(f"Prepared {len(ligand_pdbqt_files)} ligands for docking.")
    return ligand_pdbqt_files


### Protein Preparation ###

def fix_pdb_structure(input_pdb_path, output_dir):
    """Fix PDB structure(s)."""
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = []
    if input_path.is_dir():
        pdb_files = list(input_path.glob("*.pdb"))
    elif input_path.suffix.lower() == ".pdb":
        pdb_files = [input_path]
    else:
        raise ValueError(f"No PDB files found at {input_path}")
    
    fixed_paths = []
    for pdb in pdb_files:
        fixer = PDBFixer(filename=str(pdb))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(False)
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.4)

        output_pdb = output_dir / f"{pdb.stem}.pdb"
        with open(output_pdb, "w") as out:
            PDBFile.writeFile(fixer.topology, fixer.positions, out)
        
        logger.info(f"Saved fixed protein structure: {output_pdb}")
        fixed_paths.append(output_pdb)
    
    return fixed_paths


def prepare_protein(pdb_files, output_dir):
    """Convert fixed PDB(s) to PDBQT format."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(pdb_files, (str, Path)):
        pdb_files = [Path(pdb_files)]
    
    pdbqt_paths = []
    for pdb in pdb_files:
        pdbqt_file = output_dir / f"{pdb.stem}.pdbqt"
        cmd = f"obabel -ipdb {pdb} -h --partialcharge gasteiger -O {pdbqt_file} -xr"
        run(cmd)
        logger.info(f"Generated PDBQT: {pdbqt_file}")
        pdbqt_paths.append(pdbqt_file)
    
    return pdbqt_paths


### Docking ###

def run_docking(
    receptors, 
    ligands, 
    center_coords, 
    box_size, 
    smina_path="smina", 
    exhaustiveness=8, 
    num_modes=10, 
    output_dir="docked"
):
    """Dock ligand(s) into protein(s)."""
    os.makedirs(output_dir, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(receptors, (str, Path)): receptors = [receptors]
    if isinstance(ligands, (str, Path)): ligands = [ligands]

    results = []

    for rec in receptors:
        rec = Path(rec)
        for lig in ligands:
            lig = Path(lig)
            out = output_dir / f"{rec.stem}_{lig.stem}_docked.pdbqt"
            
            cmd = [
                smina_path,
                "--receptor", str(rec),
                "--ligand", str(lig),
                
  
        "--center_x", str(center_coords[0]) "
        f"--center_y {center_coords[1]} "
        f"--center_z {center_coords[2]} "
        f"--size_x {box_size[0]} "
        f"--size_y {box_size[1]} "
        f"--size_z {box_size[2]} "
        f"--num_modes, {num_modes} "
        f"--exhaustiveness {exhaustiveness} "
        f"--out {docked_pdbqt} "
    )
    
    run(cmd, description=f"Docking ligand {os.path.basename(ligand)}")
    logger.info(f"Docked ligand saved: {docked_pdbqt}")
    return docked_pdbqt

