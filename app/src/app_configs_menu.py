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
    min_samples = st.number_input("Minimum number of samples per haplotype for analysis",
                                  help = "Sometimes we get genes with a huge number of rare haplotypes, which can make interpreting plots difficult. To prevent this, plots only show data for a haplotype if the number of samples with that haplotype exceeds a threshold. We recommend a threshold of 25, but you can investigate rare haplotypes by lowering this threshold. ",
                                  min_value = 1, value = 25)
    
    toast_message = f"Threshold for minimum number of samples per haplotype for analysis has been changed to {min_samples}."

    if "min_samples" not in st.session_state:
        st.session_state["min_samples"] = min_samples

    elif st.session_state["min_samples"] != min_samples:
        st.toast(toast_message)
        st.session_state["min_samples"] = min_samples

    return min_samples

def _config_data_statistics_section(min_samples, df_haplotypes, df_join, gene_id_selected):
    _process_gene_facts(min_samples, df_haplotypes, df_join, gene_id_selected)
    return

def _config_download_data_section(gene_id_selected, df_haplotypes, df_join):

    @st.cache_data
    def _encode_df(df):
        return df.to_csv().encode('utf-8')
    
    st.markdown("Hover over the download buttons for more information on the data.")

    # Hacky quickfix to the dataframe exporting
    st.download_button("Download population-level summary",
                       _encode_df(df_haplotypes.reset_index(drop = True)),
                       file_name = f'pf-haploatlas-{gene_id_selected}_population_summary.csv',
                       help = '''Explanation of columns: "ns_changes" describes the amino acid changes of each unique haplotype; "number_of_mutations" describes number of mutations relative to 3D7; "SA", "AF-W", "AF-C", etc. shows number of samples observed with that haplotype in each geographic distribution (see sidebar for details); "Total" is the total number of samples with that haplotype; "ns_changes_list" is a list of amino acid changes of the haplotype; "sample_names" describes which lab strains the haplotype is found in''',
                       use_container_width = True)
    
    st.download_button("Download sample-level summary",
                       _encode_df(df_join.drop(columns = ["index"])),
                       file_name = f'pf-haploatlas-{gene_id_selected}_sample_summary.csv',
                       help = '''Explanation of columns: "Exclusion reason" describes the reason for a sample's removal from analysis;	"ns_changes" describes the amino acid changes of the sample for the gene selected; "Sample" is the sample name; "Study" is the clinical study of origin; "Country" of sample collection; "Admin level"	is the location of sample collection; "latitude", "longitude, "Year" of sample collection; "ENA" is the ID in the European Nucleotide Archive; "All samples same case" is reformatted sample name, "Population" refers to geographic distribution (see sidebar for details); "% callable" of SNPs, "QC pass" is whether the sample passed quality control for Pf7, "Sample type" for sequencing, "Sample was in Pf6" is whether the sample was in the previous Pf6 data resource''',
                       use_container_width = True)
    return
    
def _config_plot_settings_section():
    sample_count_mode = st.radio("Select y-axis mode for Haplotype UpSet plot:", 
                                 ["Sample counts", "Sample counts on a log scale"],
                                 index = 0)
    return sample_count_mode

def _process_gene_facts(min_samples,
                        df_haplotypes,
                        df_join,
                        gene_id_selected,
                        job_logs_file="app/files/job_logs.json"):

    with open(job_logs_file, "r") as file:
        job_logs = json.load(file)
    
    gene_info = job_logs[gene_id_selected]

    # Extract statistics
    pf7_qc_pass            = len(df_join.loc[df_join['QC pass']==True])
    missing_genotype_calls = gene_info.get('c_missing', 'N/A')
    heterozygous_calls     = gene_info.get('c_het_calls', 'N/A')
    stop_codons            = gene_info.get('c_stop_codon', 'N/A')
    sample_below_threshold = df_haplotypes.loc[df_haplotypes['Total'] < min_samples].Total.sum()
    excluded_samples       = int(missing_genotype_calls + heterozygous_calls + stop_codons + sample_below_threshold)
    included_samples       = int(pf7_qc_pass - excluded_samples)

    statistics = {
        "Sample missing genotype call":                missing_genotype_calls,
        "Heterozygous sample":                         heterozygous_calls,
        "Stop codon found for gene":                   stop_codons,
        f"Less than sample threshold ({min_samples})": sample_below_threshold
    }

    statistics_table = pd.DataFrame(list(statistics.items()),
                                    columns = ['Exclusion Reason', 'Sample Count'])
    
    statistics_table["Percentage"] = statistics_table["Sample Count"].apply(lambda x: '{:.1f} %'.format(100 * x / pf7_qc_pass))

    st.write(f"{included_samples} samples ({included_samples / pf7_qc_pass * 100:.1f} %) are available for analysis out of {pf7_qc_pass} QC pass samples for this gene, after {excluded_samples} samples ({excluded_samples / pf7_qc_pass * 100:.1f} %) have been excluded for the reasons listed below. ")
    st.dataframe(statistics_table, use_container_width = True, hide_index = True)

    return
