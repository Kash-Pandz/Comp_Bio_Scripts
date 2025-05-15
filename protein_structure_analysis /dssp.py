import mdtraj as md
import pandas as pd

def get_dssp(pdb_file):
    """Compute DSSP (secondary structure) for an input PDB file."""

    # Load the PDB file
    traj = md.load(pdb_file)
    
    # Compute DSSP secondary structure
    dssp = md.compute_dssp(traj)[0]
    
    # Extract residue information
    topology = traj.topology
    residues = [(res.resSeq, res.name) for res in topology.residues]
    
    # Create DataFrame with residue information and DSSP
    df = pd.DataFrame(residues, columns=['resid', 'resn'])
    df['dssp'] = dssp
    
    return df
