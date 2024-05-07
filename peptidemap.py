class PeptideMap():
    
    def __init__(self, protein_id, protein_file_path, peptide_file_path, strains, save_path, output_file): 
        
        protein_sequence = self.get_protein_sequence(protein_id, protein_file_path)
        strain_peptides = self.get_peptide_sequence(protein_id, peptide_file_path, strains)
        self.peptide_location(protein_id, strain_peptides, protein_sequence, output_file)
        self.plot_peptide_location(protein_id,protein_sequence,strain_peptides,strains,save_path)

    def get_protein_sequence(self, protein_id, protein_file_path):
        
        protein_id = str(protein_id)
        extract_protein = ""  
        # Set flag to indicate protein found
        found_protein = False   
    
        # Open the file and read line by line
        with open(protein_file_path, 'r') as file:
            for line in file:
                line = line.strip()                
                if line.startswith(">") and protein_id in line:                   
                    found_protein = True 
                elif found_protein:
                    if line.startswith(">") and protein_id not in line:
                        found_protein = False
                        break                    
                    elif not line.startswith(">"):
                        extract_protein += line

            cleaned_sequence = extract_protein
        return cleaned_sequence
    

    def get_peptide_sequence(self, protein_id, peptide_file_path, strains):
        # read the peptide files
        df = pd.read_excel(peptide_file_path)
       
        # create a Dataframe to store the number of peptides 
        peptide_number = pd.DataFrame()

        # create a dictionary to store the peptides from different strains
        strain_peptides={}
        
        # store the location of different strains 
        column_name = df.columns.tolist()
        
        strain_location = [column_name[3:6],column_name[6:9]]
        

        for i in range(len(strains)):
            strain_peptides[strains[i]] = [] 

            # select the rows with target protein_id
            identify_protein = df[df['Gene Name']== protein_id]
            
            
            # remove the rows with value 0
            filter_protein = identify_protein[~identify_protein[strain_location[i]].isnull().all(axis=1)]
            
            peptides = filter_protein['Sequence'].tolist()
        
            identify_protein[strains[i]] = identify_protein[strain_location[i]].sum(axis=1)
        
            strain_peptides[strains[i]].extend(peptides)
        
            peptide_number['Protein_Name'] = identify_protein['Gene Name']
            peptide_number['Sequence'] = identify_protein['Sequence']
            peptide_number[strains[i]] = identify_protein[strains[i]]
        
        output_file_pro = f'peptide_{protein_id}.csv'
        peptide_number.to_csv(output_file_pro, index=False)
             
        return strain_peptides

    def peptide_location(self, protein_id, strain_peptides, protein_sequence, output_file):
        protein_sequence = Seq(protein_sequence)
        location =[]
        for key, peptides in strain_peptides.items():
            for peptide in peptides:
                peptide = Seq(peptide)
                start_pos = protein_sequence.find(peptide)
                start_pos_n = start_pos + 1
                end_pos = start_pos_n + len(peptide) - 1
                location.append([protein_id, key, str(protein_sequence), str(peptide), start_pos_n, end_pos])
            
        output_file_loc = f'{output_file}_{protein_id}.csv'
            
        with open(output_file_loc, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Protein_Name','Strain','Protein_sequence','Peptide_sequence', 'Start_Position', 'End_Position'])
            writer.writerows(location)
    

    def plot_peptide_location(self, protein_id, protein_sequence, strain_peptides, strains, save_path = None):
        
        fig, axs = plt.subplots(len(strains), 1, figsize=(10, 10))  # Set figure size as needed
        fig.suptitle(f'{protein_id}')
        
        for i, ax in enumerate(axs):  
            ax.set_ylabel(strains[i], rotation=0)
            ax.set_yticks([])
            for peptide in strain_peptides[strains[i]]:
                start_pos = protein_sequence.find(peptide)
                start_pos = start_pos + 1
                end_pos = start_pos + len(peptide) -1
                ax.axvspan(start_pos, end_pos,ymin=-0.2, ymax=1.2, facecolor='orange', alpha=1, zorder=3)

            ax.set_xlim(0, len(protein_sequence) + 1)  # Set x-axis limits
            ax.set_ylim(-1, 1)  # Set y-axis limits
            ax.margins(y=0.1)  # Add margin to y-axis

            # Remove ticks and labels from x-axis
            ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

            # Remove outside box
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            # Display x-axis as a real line
            ax.axhline(0, color='black', linewidth=0.5)

           # Add grey box for peptide locations
            # for peptide in strain_peptides[strains[i]]:
            #     start_pos = protein_sequence.find(peptide)
            #     start_pos = start_pos +1
            #     end_pos = start_pos + len(peptide) - 1
                #ax.text((start_pos + end_pos + 1) / 2, -0.5, peptide, ha='center', va='center', color='gray', fontsize=10,
                       #bbox=dict(facecolor='lightblue', alpha=1.0, boxstyle='round,pad=0.2', linewidth=0))  # Make box opaque
        if save_path:
            plt.savefig(f"{save_path}_{protein_id}.png",dpi=300, bbox_inches='tight')  # Save plot with specified filename
        else:
            plt.show()  
            

"""import libraries"""
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import pandas as pd
import csv


