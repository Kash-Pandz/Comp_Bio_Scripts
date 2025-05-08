from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO, Select


class ChainSelect(Select):
    """Custom selector to extract specific chains"""

    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids
    

def cif_to_pdb(input_dir, chain_ids=['A']):
    """Converts cif file to pdb file in directory."""

    input_dir = Path(input_dir)
    cif_files = input_dir.glob("*.cif")

    for cif_file in cif_files:
        pdb_file = cif_file.with_suffix(".pdb")
        
        # Read .cif files
        structure = MMCIFParser.get_structure(cif_file.stem, str(cif_file))

        # Write .pdb files, selecting the specified chains
        io = PDBIO.set_structure(structure)
        io.save(str(pdb_file), select=ChainSelect(chain_ids))

if __name__ == "__main__":

    input_dir = "path/to/cif/files"
    cif_to_pdb(input_dir, chain_ids=['A'])
