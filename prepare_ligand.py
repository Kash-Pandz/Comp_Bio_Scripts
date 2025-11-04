from pathlib import Path
import os
import tempfile
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from utils import run, logger


def prepare_ligand(smiles_file, output_dir):
    """Convert ligand SMILES to PDBQT files."""
    output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    ligand_pdbqt_files = []

    with open(smiles_file) as f:
        for i, line in enumerate(f):
            fields = line.strip().split()
            if not fields:
                continue

            smi = fields[0]
            lig_name = fields[1] if len(fields) > 1 else f"ligand_{i}"
            cleaned_name = "".join(c if c.isalnum() or c in "_-" else "_" for c in lig_name)
            short_name = (cleaned_name[:3].upper()).ljust(3, "-")

            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                logger.warning(f"Skipping invalid SMILES: {smi}")
                continue
            
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)

            with tempfile.NamedTemporaryFile(suffix=".mol", delete=False) as tmp:
                Chem.MolToMolFile(mol, tmp.name)
                tmp_path = tmp.name

            pdbqt_path = output_dir / f"{short_name}.pdbqt"

            # Convert using Open Babel
            try:
                run(
                    f"obabel {tmp_path} -O {pdbqt_path} --gen3d --partialcharge gasteiger",
                    description=f"Ligand {lig_name}"
                )
            except subprocess.CalledProcessError as e:
                    logger.error(f"Open Babel failed for {lig_name}: {e}")
                    os.remove(tmp_path)
                    continue
            finally:
               # Clean up temporary file
               if os.path.exists(tmp_path):
                    os.remove(tmp_path)

            # Patch residue naem and add REMARK
            with open(pdbqt_path, "r+") as f_out:
                lines = f_out.readlines()
                new_lines = [f"REMARK Name: {short_name}\n"]
                for line in lines:
                    if line.startswith(("HETATM", "ATOM")):
                        line = line[:17] + short_name + line[20:]
                    new_lines.append(line)
                f_out.seek(0)
                f_out.writelines(new_lines)
                f_out.truncate()
            
            ligand_pdbqt_files.append(str(pdbqt_path))
                    
        
    logger.info(f"Prepared {len(ligand_pdbqt_files)} ligands for docking.")
    return ligand_pdbqt_files