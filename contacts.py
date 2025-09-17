import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import pandas as pd


def find_contacts(pdb_file, selection1, selection2, cutoff_dist):
    """Find protein-ligand or protein-protein contacts using a cutoff distance from a pdb file."""
    
    universe = mda.Universe(pdb_file)
    
    group1 = universe.select_atoms(selection1)
    group2 = universe.select_atoms(selection2)
    
    # Calculate contacts between two atom groups using a cutoff distance
    contacts = capped_distance(g1, g2, cutoff_dist, return_distances=False)
    left, right = contacts.T

    unique_res_pairs = []
    seen = set()
    for i, (l,r) in enumerate(zip(left, right)):
        segid1, resname1, resid1 = g1[l].segid, g1[l].resname, g1[l].resid
        segid2, resname2, resid2 = g2[r].segid, g2[r].resname, g2[r].resid
        dist = distances[i]

        pair_key = tuple(sorted([(segid1, resid1), (segid2, resid2)]))
        if pair_key in seen:
            continue
        seen.add(pair_key)

    unique_res_pairs.append({
            "segid1": segid1, "resid1": resid1, "resname1": resname1,
            "segid2": segid2, "resid2": resid2, "resname2": resname2,
            "distance": dist
        })
    
    
    return pd.DataFrame(data, columns=[
        "segid1", "resid1", "resname1",
        "segid2", "resid2", "resname2",
        "distance"
    ])
