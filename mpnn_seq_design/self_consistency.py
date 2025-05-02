import subprocess
import argparse
from pathlib import Path
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms


def get_tmscore(ref_pdb: Path, sel_pdb: Path, tmscore_bin: Path) -> float:
    """Calculate the TM-Score beteen two PDB files."""
    
    tmscore_cmd = subprocess.run([str(tmscore_bin), str(ref_pdb), str(sel_pdb)],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    for line in tmscore_cmd.stdout.splitlines():
        if line.startswith("TM-score"):
            # Extract the TM-score value
            tm_score_str = line.split("=")[1].strip().split()[0]
            tm_score = float(tm_score_str)
            return tm_score
    return 0.0    


def get_rmsd(ref_pdb: Path, sel_pdb: Path) -> float:
    """Calculate the RMSD between two PDB files."""
    
    # Load reference and selection pdb files into MDAnalysis universe
    ref_universe = mda.Universe(ref_pdb)
    sel_universe = mda.Universe(sel_pdb)

    # Select atom groups
    ref_group = ref_universe.select_atoms("protein and name CA")
    sel_group = sel_universe.select_atoms("protein and name CA")

    assert len(ref_group) == len(sel_group), "Structures have different number of atoms."

    # Calculate RMSD between two structures
    rmsd = rms.rmsd(ref_group.positions,
                    sel_group.positions,
                    center=True,
                    superposition=True
                   )
    
    return rmsd


def get_metrics(ref_pdb: Path, sel_dir: Path, tmscore_bin: Path, out_fname: Path):
    """Collect TM-Scores and RMSD."""
    metrics = []
    pdb_files = list(sel_dir.glob("*.pdb"))
    for sel_pdb in pdb_files:
                tmscore = get_tmscore(ref_pdb, sel_pdb, tmscore_bin)
                rmsd = get_rmsd(ref_pdb, sel_pdb)
                metrics.append({"ref_pdb": ref_pdb.stem, "sel_pdb": sel_pdb.stem, "tmscore": tmscore, "rmsd": rmsd})
    
    df = pd.DataFrame(metrics)
    df.to_csv(out_fname, index=False)


def main():
    
    parser = argparse.ArgumentParser(description="Calculate self-consistency TM-score and RMSD between reference and selection models.")
    parser.add_argument('--ref_pdb', required=True, help="Path to the reference PDB.")
    parser.add_argument('--sel_dir', required=True, help="Directory with selection PDB files.")
    parser.add_argument('--tmscore_bin', required=True, help="Path to the TM-score binary.")
    parser.add_argument('--out_fname', required=True, help="Output file .csv format.")

    args = parser.parse_args()

    get_metrics(Path(args.ref_pdb,), Path(args.sel_dir), Path(args.tmscore_bin), Path(args.out_fname))

if __name__ == "__main__":
      main()

