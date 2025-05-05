import os
import glob
import argparse
from tqdm import tqdm
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms


ATOM_SELECTIONS = {
    'calpha': 'name CA',
    'backbone': 'backbone',
    'heavy': 'not name H*'
}

def select_atom_group (universe, atom_type, chain_id=None, custom=None):
    """Create an atom group selection string."""
    base_sel = ATOM_SELECTIONS.get(atom_type)
    if not base_sel:
            raise ValueError(f"Invalid atom atype: {atom_type}")
    
    # Combine base selection if custom selection provided
    selection = f"{base_sel} and ({custom})" if custom else base_sel

    if chain_id:
        selection += f" and segid {chain_id}"
    
    return universe.select_atoms(selection)


def get_rmsd(ref_atoms, sel_atoms):
    """Calculate pairwise rmsd between two protein structures."""
    if len(ref_atoms) != len(sel_atoms):
        raise ValueError("Atom count mismatch between ref_atoms and sel_atoms.")
    return rms.rmsd(ref_atoms.positions, sel_atoms.positions, center=True, superposition=True)


def main(ref_pdb, pdb_dir, atom_type, chain_id, custom_selection, output_csv):
    ref_universe = mda.Universe(ref_pdb)
    pdb_files = glob.glob(os.path.join(pdb_dir, '*.pdb'))

    if not pdb_files:
        print("No PDB files found.")
        return

    ref_atoms = select_atom_group(ref_universe, atom_type, chain_id, custom_selection)
    results = []

    for pdb_path in tqdm(pdb_files, desc="Processing PDBs"):
        pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
        target_universe = mda.Universe(pdb_path)
        target_atoms = select_atom_group(target_universe, atom_type, chain_id, custom_selection)
        rmsd_val = get_rmsd(ref_atoms, target_atoms)
        results.append({'reference': os.path.basename(ref_pdb), 'query': pdb_name, 'rmsd': rmsd_val})

    pd.DataFrame(results).to_csv(output_csv, index=False)
    print(f"RMSD results saved to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pairwise RMSD calculator for PDB files.")
    parser.add_argument('--ref_pdb', required=True, help='Reference PDB file')
    parser.add_argument('--pdb_dir', required=True, help='Directory of PDB files')
    parser.add_argument('--atom_type', choices=['calpha', 'backbone', 'heavy'], default='calpha')
    parser.add_argument('--chain_id', help='Optional chain id')
    parser.add_argument('--custom_selection', help='Additional selection string (e.g., resid 45:50)')
    parser.add_argument('--output_csv', default='rmsd_results.csv')

    args = parser.parse_args()
    main(args.ref_pdb, args.pdb_dir, args.atom_type, args.chain_id, args.custom_selection, args.output_csv)
