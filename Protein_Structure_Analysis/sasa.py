import mdtraj as md
import pandas as pd


## Max SASA values from Tien et al. 2013 ##

MAX_SASA= {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0,
    'GLN': 223.0, 'GLU': 225.0, 'GLY': 104.0, 'HIS': 224.0, 'ILE': 197.0,
    'LEU': 201.0, 'LYS': 236.0, 'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0,
    'SER': 155.0, 'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0
}


def get_sasa(pdb_file):
    """Calculate the sasa and the rsa of each residue in pdb file."""
    traj = md.load(pdb_file)
    sasa = md.shrake_rupley(traj) * 100
    topology = traj.topology

    data = []

    for residue in topology.residues:
        atom_indices = [atom.index for atom in residue.atoms]
        sasa_value = np.sum(sasa[0, atom_indices])  
        
        res_name = residue.name

        # Get max SASA
        max_sasa = MAX_SASA.get(res_name, None)
        
        # Calculate the relative SASA for each residue
        rsa_value = sasa_value / max_sasa if max_sasa else None
        
        data.append([residue.index + 1, res_name, sasa_value, rsa_value])

    # Store results in pandas dataframe
    df = pd.DataFrame(data, columns=["resid", "resname", "sasa ($\AA$)", "rsa"])
    return df
