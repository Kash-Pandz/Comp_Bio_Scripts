import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
import matplotlib.pyplot as plt



def plot_bb_dihedrals(pdb_file):
    """Plot the phi and psi backbone dihedral angles from a PDB file."""

    universe = mda.Universe(pdb_file)
    protein = universe.select_atoms("protein")
    rama = dihedrals.Ramachandran(protein).run()

    fig = rama.plot(color='black', marker='.', ref=True)
    plt.show()



def plot_sc_dihedrals(pdb_file):
    """Plot the X1 and X2 dihedral angles from a PDB file."""

    universe = mda.Universe(pdb_file)
    protein = universe.select_atoms("protein")
    sc = dihedrals.Janin(protein).run()

    fig = sc.plot(color='black', marker='.', ref=True)
    plt.show()



  
