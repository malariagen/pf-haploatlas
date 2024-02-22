import streamlit as st
import json
from src.utils import cache_load_population_colours, _cache_load_utility_mappers

def process_gene_facts(gene_id_selected, job_logs_file="path/to/job_logs.json"):
    """Main function called in main.py to show gene facts"""
    
    utility_mappers = _cache_load_utility_mappers()

    # Load job logs from the JSON file
    with open(job_logs_file, "r") as file:
        job_logs = json.load(file)
    
    # Initialize gene_info with an empty dictionary
    gene_info = {}

    # Check if gene_id_selected is present in job_logs
    if gene_id_selected in job_logs:
        gene_info = job_logs[gene_id_selected]
        
    # Display gene facts using job_logs information
    st.subheader(f'Viewing gene: {utility_mappers["gene_ids_to_gene_names"][gene_id_selected]}')

    st.write(f"While {gene_info.get('c_inc_s', 'N/A')} ({gene_info.get('c_inc_s', 'N/A')/16103*100:.2f}%) QC-passed Pf7 samples included in the analysis,")
    st.write(f"{gene_info.get('c_missing', 'N/A')} ({gene_info.get('c_missing', 'N/A')/16103*100:.2f}%) samples with missing genotype calls,")
    st.write(f"{gene_info.get('c_het_calls', 'N/A')} ({gene_info.get('c_het_calls', 'N/A')/16103*100:.2f}%) samples with heterozygous calls, and")
    st.write(f"{gene_info.get('c_stop_codon', 'N/A')}{gene_info.get('c_stop_codon', 'N/A')/16103*100:.2f}%) samples with stop codons excluded from the analysis.")