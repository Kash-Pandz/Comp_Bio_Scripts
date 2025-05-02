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


def get_pae(pdb_dir):
    """Obtain the Predicted Alignment Error (PAE) for AF predictions."""

    pae = []

    files = sorted(glob.glob(os.path.join(pdb_dir, "*.json")))

    for file in files:
        with open(file, 'r') as f:
            json_data = json.load(f)
            pae_values = json_data.get('pae', [])
            if isinstance(pae_values, list) and pae_values:
                mean_pae = np.mean(pae_values)
                pae.append({"pae": mean_pae})
    
    return pae


def get_AF_metrics(plddt, pae):
    """Compile the plddt and pae metrics into a pandas dataframe."""
    
    assert len(plddt) == len(pae), "Metrics list must be the same length"

    data = [
        {**plddt_data, **pae_data} for plddt_data, pae_data in zip(plddt, pae)
    ]

    return pd.DataFrame(data)
