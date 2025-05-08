import pandas as pd
import MDAnalysis as mda


def pdb_to_df(pdb_file: str) -> pd.DataFrame:
    """Convert PDB file into pandas DataFrame."""
    
    universe = mda.Universe(pdb_file)

    data = []
    
    for atom in universe.atoms:
        coords = atom.position
        data.append({
            "chain": atom.segid,
            "resname": atom.resname,
            "res_id": atom.resid,
            "atom_id": atom.index,
            "atom_name": atom.name,
            "x": round(coords[0], 3),
            "y": round(coords[1], 3),  
            "z": round(coords[2], 3),  
            "bfactor": atom.tempfactor if atom.tempfactor is not None else None 
        })
    
    df = pd.DataFrame(data)
    
    return df
