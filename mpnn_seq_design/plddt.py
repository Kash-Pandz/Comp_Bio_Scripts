import os
import glob
import json
import numpy as np
import MDAnalysis as mda
import pandas as pd



def get_plddt(pdb_dir):
    """Obtain the predicted local distance difference test (pLDDT) for AF predictions."""

    plddt = []

    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))

    for pdb_file in pdb_files:
        universe = mda.Universe(pdb_file)
        mean_plddt = np.mean(universe.atoms.tempfactors)
        plddt.append({
            "pdb_fname": os.path.basename(pdb_file),
            "plddt": mean_plddt
        })

    return plddt
