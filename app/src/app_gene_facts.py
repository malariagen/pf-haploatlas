import streamlit as st
import json
from src.utils import cache_load_population_colours, _cache_load_utility_mappers
import pandas as pd

def process_gene_facts(min_samples, df_haplotypes, gene_id_selected, job_logs_file="path/to/job_logs.json"):
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
        qc_fail = 4661
        missing_genotype_calls = gene_info.get('c_missing', 'N/A')
        heterozygous_calls = gene_info.get('c_het_calls', 'N/A')
        stop_codons = gene_info.get('c_stop_codon', 'N/A')
        sample_below_threshold = (df_haplotypes['Total'] < min_samples).sum()
        excluded_samples = 20864 - included_samples + sample_below_threshold
        # Create data for the table
        table_data = [
            {"Exclusion Reason": "QC fail", "Count": qc_fail, "Percentage": qc_fail / 20864 * 100},
            {"Exclusion Reason": "Missing genotype calls", "Count": missing_genotype_calls, "Percentage": missing_genotype_calls / 20864 * 100},
            {"Exclusion Reason": "Heterozygous calls", "Count": heterozygous_calls, "Percentage": heterozygous_calls / 20864 * 100},
            {"Exclusion Reason": "Stop codons", "Count": stop_codons, "Percentage": stop_codons / 20864 * 100},
            {"Exclusion Reason": f"Sample threshold ({min_samples})", "Count": sample_below_threshold, "Percentage": sample_below_threshold / 20864 * 100},
            {"Exclusion Reason": "Total", "Count": excluded_samples, "Percentage": excluded_samples / 20864 * 100}
        ]
        df_table = pd.DataFrame.from_dict(table_data)
        df_table['Percentage'] = df_table['Percentage'].apply(lambda x: '{:.2f}%'.format(x))
        print(df_table)
        st.table(df_table.set_index('Exclusion Reason', drop=True))