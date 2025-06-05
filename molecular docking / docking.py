import subprocess
import MDAnalysis as mda



def split_pdb(pdb_file, ligand_name, protein_file="protein.pdb", ligand_file="ligand.pdb"):
    """Splits input PDB file into protein and ligand files."""

    # Load PDB file into MDAnalysis universe
    universe = mda.Universe(pdb_file)

    # Select protein and ligand atoms
    protein = universe.select_atoms("protein")
    ligand = universe.select_atoms(f"resname {ligand_name}")

    # Write protein and ligand pdb files
    protein.write(protein_file)
    ligand.write(ligand_file)
    

### fpocket command ###
## fpocket -f PDB_FILE ##

def extract_center(pdb_file, selection, output_file="centre.txt"):
    """Extracts the geometric center of binding pocket or native ligand from a PDB file."""
    
    # Load PDB file into MDAnalyisis universe
    universe = mda.Universe(pdb_file)

    # Select ligand atoms 
    ligand = universe.select_atoms(f"resname {ligand_name}")

    if len(ligand) > 0:
        cog = ligand.center_of_geometry()
        name = f"Ligand ({selection})"
    else:
        pocket = universe.select_atoms("all")  # Selection made after manual inspection of pocket PDB file obtained from fpocket
        if len(pocket) > 0:
            cog = pocket.center_of_geometry()
            name = "pocket"
        else:
            print("No ligand or pocket found in the PDB file.")
            return None

     return cog.tolist()


def protein_to_pdbqt(prot_pdb, prot_pdbqt):
   """Convert protein pdb into a pdbqt file."""
    command = [
        'obabel', 
        '-ipdb', prot_pdb,
        '-h',
        '--partialcharge', 'eem',
        '-opdbqt', 
        '-O', prot_pdbqt, 
        '-xr'
    ]
    subprocess.run(command, stdout=subprocess.DEVNULL)


def ligand_to_pdbqt(ligand_pdb, ligand_pdbqt):
   """Convert ligand sdf or pdb file into pdbqt format."""
    command = [
        'obabel', 
        '-ipdb', lig_pdb,
        '-h',
        '--partialcharge', 'eem',
        '-opdbqt', 
        '-O', lig_pdbqt, 
    ]
    subprocess.run(command, stdout=subprocess.DEVNULL)


def run_docking(protein_pdbqt, ligand_pdbqt, center_coords, box_size, poses, exhaustiveness=8, docked_pdbqt="docked_poses.pdbqt"):
    """Run SMINA molecular docking."""

    smina_cmd = [
        "smina",
        "--receptor", protein_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center_coords[0]),
        "--center_y", str(center_coords[1]),
        "--center_z", str(center_coords[2]),
        "--size_x", str(box_size[0]), 
        "--size_y", str(box_size[1]), 
        "--size_z", str(box_size[2]),
        "--num_modes", str(poses),
        "--exhaustiveness", str(exhaustiveness), 
        "--out", docked_pdbqt,
    ]
    result = subprocess.run(smina_cmd, capture_output=True, text=True)
    

