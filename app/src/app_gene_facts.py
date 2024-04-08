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
    # little trick to adjust the width of expander, since markdown css can not be specifed just for the second expander
    col1, col2 = st.columns([3, 2])
    with col1.expander("Click to view how many Pf7 samples are included/excluded."):
        included_samples = gene_info.get('c_inc_s', 'N/A')
        excluded_samples = 20864 - included_samples if included_samples != 'N/A' else 'N/A'
        qc_fail = 4661
        missing_genotype_calls = gene_info.get('c_missing', 'N/A')
        heterozygous_calls = gene_info.get('c_het_calls', 'N/A')
        stop_codons = gene_info.get('c_stop_codon', 'N/A')
        st.write(f"**Included:** {included_samples} ({included_samples/20864*100:.2f}%)")
        st.write(f"""
            **Excluded:** {excluded_samples} ({excluded_samples/20864*100:.2f}%)
            | Exclusion Reason           | Count      | Percentage |
            |----------------------------|------------|------------|
            | QC fail                    | {qc_fail}   | 22.34%     |
            | Missing genotype calls     | {missing_genotype_calls} | {missing_genotype_calls/20864*100:.2f}% |
            | Heterozygous calls         | {heterozygous_calls} | {heterozygous_calls/20864*100:.2f}% |
            | Stop codons                | {stop_codons} | {stop_codons/20864*100:.2f}% |
            | Total                 | {excluded_samples} | {excluded_samples/20864*100:.2f}% |
        """)