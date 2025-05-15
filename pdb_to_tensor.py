import MDAnalysis as mda
import torch


def convert_pdb_to_tensor(pdb_file):
    """ A funtion to convert PDB atom coordinates to pytorch tensor."""

    universe = mda.Universe(pdb_file)

    # Select C-alpha atoms and get their positions
    ca_atoms = u.select_atoms("name CA")
    ca_tensor = torch.tensor(ca_atoms.positions, dtype=torch.float32)
    print(ca_tensor.size())
    
    # Select backbone atoms (N, CA, C, O) and get their positions
    backbone_atoms = u.select_atoms("backbone")
    backbone_tensor = torch.tensor(backbone_atoms.positions, dtype=torch.float32)
    print(backbone_tensor.size())

    return ca_tensor, backbone_tensor
