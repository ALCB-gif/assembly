#Import libraries
from libdl.dl import dl
import pandas as pd 
import os
import numpy as np
import pandas as pd 
from libgreen.models import Genome
from libgreen.models import Assembly
from libgreen.models import Protein

###### Non spores ######
"""
How to use me:

Run non_spores.py followed by 
1. Location for assemblies
2. Location for proteins 

If a location is not given it will run in the place you
are running the code from...
"""
#Start sql query
gcc_data = pd.read_sql("""
SELECT * FROM ds_mr_global.gcc_select
""", con=dl.engine)

#SQL query (risk group)
gcc_data = pd.read_sql("""
SELECT * FROM ds_mr_biohealth.gspo_strains where
qps_amr = 'The strains should not harbour any acquired antimicrobial resistance genes to clinically relevant antimicrobials.'
and qps_toxin = 'Absence of toxigenic activity.'
and qps_production_only is NULL
and riskgroup= '1' 
""", con=dl.engine)

gcc_data = pd.read_sql("""
SELECT * FROM ds_mr_biohealth.lab_gcc_biohealth_insilico_summary
""",con=dl.engine)

#Data preview 
gcc_data.to_csv('all_data_s.tsv', sep='\t', index=False)

#Output directories
g_output_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.getcwd(), "assembly")
p_output_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.getcwd(), "proteins")

#Create directories if non existent
os.makedirs(g_output_path, exist_ok=True)
os.makedirs(p_output_path, exist_ok=True)

protein_list = {}
for asm_id in gcc_data["assembly_identifier"]:
    if pd.notnull(asm_id): 
        try:
            temp_genome = Assembly.find(asm_id)
            if temp_genome:
                asm = Assembly.get(temp_genome[0].accession_number)
                
                #Strains (Genomes)
                genome_output_name = asm.organism_name.replace(" ", "_") + "_" + asm.accession_number + ".fna"
                genome_output_path = os.path.join(g_output_path, genome_output_name)
                asm.write_dna('fasta', genome_output_path)
                
                #Proteins
                proteins = Protein.get([protein.key for protein in asm.proteins()])
                protein_output_name = asm.organism_name.replace(" ", "_") + "_" + asm.accession_number + ".faa"
                protein_output_path = os.path.join(p_output_path, protein_output_name)
                
                with open(protein_output_path, "w") as f:
                    for protein in proteins:
                        f.write(f">{protein.key}\n{protein.sequence}\n")
                
                print(f"Processed assembly {asm_id}:")
                print(f" Genome file: {genome_output_path}")
                print(f" Protein file: {protein_output_path}")
        
        except Exception as e:
            print(f"Error processing asm_id {asm_id}: {e}")
