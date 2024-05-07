import os 
import peptidemap
from peptidemap import PeptideMap

current_dir = os.path.dirname(os.path.abspath(__file__))
print(current_dir)

if __name__ == "__main__":
    # target protein and strains    
    protein = ['enzyme']
    strains = ['strain1','strain2']

    protein_file_path = current_dir + "/protein_sequence.fasta"
    peptide_file_path = current_dir + "/example.xlsx"

    save_path= "peptide_map"
    output_file = "peptide_location"
    
    for protein_id in protein:
        Map = PeptideMap(protein_id, protein_file_path, peptide_file_path, strains, save_path, output_file)
    
