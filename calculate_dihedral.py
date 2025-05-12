import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran

def get_bb_dihedrals(pdb_file):
    """Calculates the phi and psi backbone dihedral angles from a PDB file."""

    universe = mda.Universe(pdb_file)
    protein = universe.select_atoms("protein")

    rama = Ramachandran(protein).run()
    angles = rama.angles[0]

    data = []
    for res, (phi, psi) in zip(rama.residues, angles):
        data.append({
            'Residue': res.resid,
            'Residue Name': res.resname,
            'Phi (°)': np.degrees(phi),
            'Psi (°)': np.degrees(psi)
        })

    # Convert the list of dictionaries to a pandas DataFrame
    df = pd.DataFrame(data)
    return df
  
  
