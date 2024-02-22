####################################
### The script identifies genes located within the core regions of the 3D7 genome.
### Iterating over the 5568 genes in 3D7:
###     For each gene, it examines the start-end coordinates of its coding sequences (CDS):
###      - Excludes the gene if none of its CDS falls within the core region
###      - Includes the gene if at least one CDS is within the core, marking it as a boundary
###      - Includes the gene if all of its CDS are entirely within the core region
####################################


import allel
import pandas as pd
import json

# Load and extract 3D7 gene annotations into a list
gff_fn = '/lustre/scratch124/gsu/legacy/pfalciparum/resources/snpEff/data/Pfalciparum_GeneDB_Feb2020/Pfalciparum_replace_Pf3D7_MIT_v3_with_Pf_M76611.gff'
df_gff = allel.gff3_to_dataframe(gff_fn, attributes=['ID', 'Name'])
filtered_df_gff = df_gff[df_gff['type'] == 'gene']
Pf3D7_genes_list = filtered_df_gff['ID'].tolist()

# Extract coding regions of 3D7 (CDS excludes UTRs)
cds_df_gff = df_gff[df_gff['type'] == 'CDS']

# Load and extract core region coordinates into list
# Core region coordinates are extracted from: ftp://ngs.sanger.ac.uk/production/malaria/pf-crosses/1.0//regions-20130225.onebased.txt
regions = pd.read_csv('../../files/regions-20130225.onebased.txt', sep='\t', names=['chrom', 'start', 'end', 'region'])
core_regions = regions[regions['region'] == 'Core']

Pf3D7_core_genes = dict()

for gene_ID in Pf3D7_genes_list:
    # Filter annotations specific to this gene 
    filtered_gene_df =  df_gff.loc[( df_gff['ID'].str.startswith(gene_ID) )  & ( df_gff['type'] == 'CDS')]
    chrom = filtered_gene_df.iloc[0]['seqid']
    gene_name = filtered_df_gff[filtered_df_gff.ID == gene_ID].iloc[0]['Name']
    
    # Filter core region specific to this chromosome 
    filtered_core = core_regions[core_regions.chrom == chrom]
    
    # Compare each gene CDS start-end coordinates and core regions 
    # Keep genes that are only in the core regions
    if not filtered_core.empty:      
        for start, end in zip( filtered_gene_df['start'], filtered_gene_df['end']):
            for index, row in filtered_core.iterrows():
                core_start = row['start']
                core_end = row['end']

                if (core_start <= start < core_end) and (end <= core_end):
                    Pf3D7_core_genes[gene_ID] = gene_name

                elif (start < core_start) and (core_start < end < core_end):
                    Pf3D7_core_genes[gene_ID] = gene_name
                    print(f"Haplotypes are likely to be partial for {gene_name}.")

                elif (core_start <= start < core_end) and (core_end < end):
                    Pf3D7_core_genes[gene_ID] = gene_name
                    print(f"Haplotypes are likely to be partial for {gene_name}.")

count_core_genes = len(Pf3D7_core_genes)
total_genes_number = len(Pf3D7_genes_list)

# Logging 
print(f'Total number of non-core region genes removed:  {total_genes_number- count_core_genes}')
print(f'Total number of core region genes:  {count_core_genes}')

# Write the dictionary to a JSON file
with open('core_genes.json', 'w') as json_file:
    json.dump(Pf3D7_core_genes, json_file, indent=2)