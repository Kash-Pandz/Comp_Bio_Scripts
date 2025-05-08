import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import pandas as pd


def find_contacts(pdb_file, selection1, selection2, cutoff_dist):
    """Find protein-ligand or protein-protein contacts using a cutoff distance from a pdb file."""
    
    universe = mda.Universe(pdb_file)
    
    group1 = universe.select_atoms(selection1)
    group2 = universe.select_atoms(selection2)
    
    # Calculate contacts between two atom groups using a cutoff distance
    contacts = capped_distance(group1,
                               group2,
                               cutoff_dist,
                               return_distances=False)
    
    left, right = contacts.T
    
    group1_segids = group1[left].segids
    group2_segids = group2[right].segids
    group1_resnames = group1[left].resnames
    group2_resnames = group2[right].resnames
    group1_resids = group1[left].resids
    group2_resids = group2[right].resids

    res_pairs = [
        {
            'segid1': segid1, 'resname1': resname1, 'resid1': resid1,
            'segid2': segid2, 'resname2': resname2, 'resid2': resid2
        }
        for segid1, resname1, resid1, segid2, resname2, resid2
        in zip(group1_segids, group1_resnames, group1_resids,
               group2_segids, group2_resnames, group2_resids)
    ]

    unique_res_pairs = []
    seen = set()
    for pair in res_pairs:
        sorted_resids = tuple(sorted([pair['resid1'], pair['resid2']]))
        if sorted_resids not in seen:
            seen.add(sorted_resids)
            unique_res_pairs.append(pair)

    df = pd.DataFrame(unique_res_pairs, columns=['resid1', 
                                                 'resname1', 
                                                 'segid1', 
                                                 'resid2', 
                                                 'resname2', 
                                                 'segid2'])
    
    return df
