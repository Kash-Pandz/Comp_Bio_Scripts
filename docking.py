import os
from utils import run, logger

def run_docking(receptor, ligand, center_coords, box_size, smina_path="smina", exhaustiveness=8, num_modes=10, output_dir="."):
    """Run molecular docking."""
    os.makedirs(output_dir, exist_ok=True)
    docked_pdbqt = os.path.join(output_dir, f"docked_{os.path.basename(ligand)}")
   
    cmd = (
        f"{smina_path} "
        f"--receptor {receptor} "
        f"--ligand {ligand} "
        f"--center_x {center_coords[0]} "
        f"--center_y {center_coords[1]} "
        f"--center_z {center_coords[2]} "
        f"--size_x {box_size[0]} "
        f"--size_y {box_size[1]} "
        f"--size_z {box_size[2]} "
        f"--num_modes, {num_modes} "
        f"--exhaustiveness {exhaustiveness} "
        f"--out {docked_pdbqt} "
    )
    
    run(cmd, description=f"Docking ligand {os.path.basename(ligand)}")
    logger.info(f"Docked ligand saved: {docked_pdbqt}")
    return docked_pdbqt