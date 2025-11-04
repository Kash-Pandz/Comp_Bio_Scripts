from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from utils import run, logger


def fix_pdb_structure(input_pdb: Path, output_pdb: Path):
    fixer = PDBFixer(filename=str(input_pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    with open(output_pdb, "w") as out:
        PDBFile.writeFile(fixer.toplogy, fixer.positions, out)
    logger.info(f"Saved fixed protein structure: {output_pdb}")
    return output_pdb


def prepare_protein(pdb_file, pdbqt_file):
    cmd = f"obabel -ipdb {pdb_file} -h --partialcharge gasteiger -O {pdbqt_file} -xr"
    run(cmd, description="Protein preparation")
    return pdbqt_file
