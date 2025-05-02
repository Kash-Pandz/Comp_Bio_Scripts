import os
import glob
import subprocess
import argparse
from loguru import logger
from Bio import SeqIO
from pathlib import Path
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="LigandMPNN Sequence Design Script")
    parser.add_argument("--pdb_dir", type=str, required=True, help="Directory containing PDB files.")
    parser.add_argument("--temps", type=float, nargs='+', required=True, help="Sampling temperatures (e.g., 0.1 0.2 0.3).")
    parser.add_argument("--fixed_residues", type=str, required=True, help="Fixed residues (e.g., 'A:15,16,17').")
    parser.add_argument("--omit_AA", type=str, default=None, help="Residues to be omitted from designs (e.g., 'C').")
    parser.add_argument("--checkpoint_ligandmpnn", type=str, default="/PATH/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt", help="Path to LigandMPNN model checkpoint.")
    parser.add_argument("--checkpoint_sc", type=str, default="/PATH/LigandMPNN/model_params/ligandmpnn_sc_v_32_002_16.pt", help="Path to LigandMPNN side chain model checkpoint.")
    parser.add_argument("--number_of_batches", type=int, default=2, help="Number of batches to run (default: 2).")
    parser.add_argument("--batch_size", type=int, default=10, help="Batch size for sequence design (default: 10).")
    return parser.parse_args()

def run_seq_design(pdb_dir, temps, fixed_residues, omit_aa, checkpoint_ligandmpnn, checkpoint_sc, number_of_batches, batch_size):
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))

    if not pdb_files:
        logger.error("No PDB files found in the specified directory!")
        return

    for pdb_file in pdb_files:
        pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]

        for temp in temps:
            out_dir = Path(f"{pdb_name}_temp_{temp}")
            out_dir.mkdir(parents=True, exist_ok=True)

            mpnn_cmd = [
                "python", "run.py",
                "--seed", "111",
                "--pdb_path", pdb_file,
                "--out_folder", str(out_dir),
                "--model_type", "ligand_mpnn",
                "--number_of_batches", str(number_of_batches),
                "--batch_size", str(batch_size),
                "--temperature", str(temp),
                "--fixed_residues", fixed_residues,
            ]

            if omit_aa:
                mpnn_cmd.extend(["--omit_AA", omit_aa])

            mpnn_cmd.extend([
                "--checkpoint_ligand_mpnn", checkpoint_ligandmpnn,
                "--checkpoint_path_sc", checkpoint_sc,
                "--pack_side_chains", "1",
                "--number_of_packs_per_design", "1",
                "--pack_with_ligand_context", "1",
                "--repack_everything", "1",
                "--ligand_mpnn_cutoff_for_score", "6.0",
            ])

            logger.info(f"Running ligand_mpnn for {pdb_name} at temperature {temp}")
            logger.info("Command: " + " ".join(mpnn_cmd))

            subprocess.run(mpnn_cmd, check=True)
            logger.info(f"ligand_mpnn finished for {pdb_name} at temperature {temp}")

            # Process FASTA files
            seqs_dir = out_dir / "seqs"
            logger.info(f"Checking directory for fasta files: {seqs_dir}")

            if seqs_dir.exists():
                logger.info(f"Listing files in {seqs_dir}: {os.listdir(seqs_dir)}")

            fasta_files = list(seqs_dir.glob("*.fa"))
            if not fasta_files:
                logger.warning(f"No fasta files found in {seqs_dir}")
                continue

            input_fasta = fasta_files[0]
            temp_str = str(temp).replace('.', '')
            output_fasta = out_dir / f"{pdb_name}_seqs_{temp_str}.fa"

            records = list(SeqIO.parse(input_fasta, "fasta"))
            with open(output_fasta, "w") as out_fname:
                for i, record in enumerate(records[1:], 1):
                    new_header = f"{pdb_name}_temp_{temp_str}_seq_{i}"
                    record.id = new_header
                    record.description = ""
                    SeqIO.write(record, out_fname, "fasta")

            logger.info(f"Fasta file {input_fasta.name} headers renamed and saved to {output_fasta.name}")

    logger.info("Sequence design and FASTA processing completed!")

if __name__ == "__main__":
    args = parse_arguments()

    logger.remove()
    logger.add(sys.stdout, level="INFO", format="{time} {level} {message}")

    logger.info("Starting sequence design and FASTA processing.")
    run_seq_design(
        args.pdb_dir,
        args.temps,
        args.fixed_residues,
        args.omit_AA,
        args.checkpoint_ligandmpnn,
        args.checkpoint_sc,
        args.number_of_batches,
        args.batch_size
    )
