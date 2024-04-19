import streamlit as st
import pandas as pd
import json

from src.utils import _cache_load_utility_mappers

def process_configs_menu(gene_id_selected, df_haplotypes, df_join):
    """Main function called in main.py to handle user config settings in the expander"""
    
    col1, col2 = st.columns([2, 1]) # column width ratios

    # ----- COL1 ----- #
    with col1.expander("Click to see more about the data"):
        st.subheader("Data filtering settings")
        min_samples = _config_data_filtering_section()
        st.divider()

        st.subheader("Data filtering statistics")
        _config_data_statistics_section(min_samples, df_haplotypes, df_join, gene_id_selected)
        st.divider()

        st.subheader("Download data")
        _config_download_data_section(gene_id_selected, df_haplotypes, df_join)
        st.divider()
        
        st.subheader("Plot settings")
        sample_count_mode = _config_plot_settings_section()

    # ----- COL2 ----- #
    selected_gene_plasmodb_url = f'https://plasmodb.org/plasmo/app/record/gene/{gene_id_selected}'
             
    col2.markdown(
        f'''<a href="{selected_gene_plasmodb_url}" style="display: inline-block;
            padding: 11px 20px; background-color: #b00023;
            color: white;
            text-align: center;
            text-decoration: none;
            font-size: 16px; border-radius:
        4px; width: 100%;">Browse Gene on PlasmoDB</a>''',
        unsafe_allow_html=True
    )
    
    return min_samples, sample_count_mode

def _config_data_filtering_section():
    min_samples = st.number_input("Minimum number of samples per haplotype for analysis", min_value = 1, value = 25)
    return min_samples

def _config_data_statistics_section(min_samples, df_haplotypes, df_join, gene_id_selected):
    _process_gene_facts(min_samples, df_haplotypes, df_join, gene_id_selected)
    return

def _config_download_data_section(gene_id_selected, df_haplotypes, df_join):

    @st.cache_data
    def _encode_df(df):
        return df.to_csv().encode('utf-8')

    st.download_button("Download population-level summary",
                       _encode_df(df_haplotypes),
                       file_name = f'{gene_id_selected}_summary.csv',
                       use_container_width = True)
    
    st.download_button("Download sample-level summary",
                       _encode_df(df_join),
                       file_name = f'{gene_id_selected}_summary.csv',
                       use_container_width = True)
    return
    
def _config_plot_settings_section():
    sample_count_mode = st.radio("Select x-axis mode for Plot 1:", 
                                 ["Sample counts", "Sample counts on a log scale"],
                                 index = 0)
    return sample_count_mode

def _process_gene_facts(min_samples,
                        df_haplotypes,
                        df_join,
                        gene_id_selected,
                        job_logs_file="work/backend/job_logs.json"):

    with open(job_logs_file, "r") as file:
        job_logs = json.load(file)
    
    gene_info = job_logs[gene_id_selected]

    # Extract statistics
    pf7_total_samples      = len(df_join)
    qc_fail                = len(df_join.loc[df_join['QC pass']==False])
    missing_genotype_calls = gene_info.get('c_missing', 'N/A')
    heterozygous_calls     = gene_info.get('c_het_calls', 'N/A')
    stop_codons            = gene_info.get('c_stop_codon', 'N/A')
    sample_below_threshold = df_haplotypes.loc[df_haplotypes['Total'] < min_samples].Total.sum()
    excluded_samples       = int(qc_fail + missing_genotype_calls + heterozygous_calls + stop_codons + sample_below_threshold)
    included_samples       = int(pf7_total_samples - excluded_samples)

    statistics = {
        "Sample failed QC":                            qc_fail,
        "Sample missing genotype call":                missing_genotype_calls,
        "Heterozygous sample":                         heterozygous_calls,
        "Stop codon found for gene":                   stop_codons,
        f"Less than sample threshold ({min_samples})": sample_below_threshold
    }

    statistics_table = pd.DataFrame(list(statistics.items()),
                                    columns = ['Exclusion Reason', 'Sample Count'])
    
    statistics_table["Percentage"] = statistics_table["Sample Count"].apply(lambda x: '{:.1f} %'.format(100 * x / pf7_total_samples))

    st.write(f"{included_samples} samples ({included_samples / pf7_total_samples * 100:.1f} %) are available for analysis out of {pf7_total_samples} samples, after {excluded_samples} samples ({excluded_samples / pf7_total_samples * 100:.1f} %) have been excluded for the reasons listed below. ")
    st.dataframe(statistics_table, use_container_width = True, hide_index = True)

    return