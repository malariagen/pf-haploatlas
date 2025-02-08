import streamlit as st
import pandas as pd
import json

from src.utils import _cache_load_utility_mappers

DATAPACK_BASE_PATH = "app/files/datapack-2025-02-08"

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
    min_samples = st.number_input(
        "Minimum number of samples per haplotype for analysis",
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
    job_logs_file = f"{DATAPACK_BASE_PATH}/auxiliary/gene_log.tsv"
    gene_log = pd.read_csv(job_logs_file, sep = "\t")
    sample_exclusion_count_columns = ["c_exc_s", "c_inc_s", "c_missing", "c_het_calls", "c_stop_codon", "c_unq_h"]
    gene_info = {col: int(gene_log.loc[gene_log.gene_id == gene_id_selected, col].values[0]) for col in sample_exclusion_count_columns}

    # Extract statistics
    n_qc_pass_samples      = len(df_join.loc[df_join['QC pass']==True])
    missing_genotype_calls = gene_info.get('c_missing', 'N/A')
    heterozygous_calls     = gene_info.get('c_het_calls', 'N/A')
    stop_codons            = gene_info.get('c_stop_codon', 'N/A')
    sample_below_threshold = df_haplotypes.loc[df_haplotypes['Total'] < min_samples].Total.sum()
    excluded_samples       = int(missing_genotype_calls + heterozygous_calls + stop_codons + sample_below_threshold)
    included_samples       = int(n_qc_pass_samples - excluded_samples)

    statistics = {
        "Sample missing genotype call":                missing_genotype_calls,
        "Heterozygous sample":                         heterozygous_calls,
        "Stop codon found for gene":                   stop_codons,
        f"Less than sample threshold ({min_samples})": sample_below_threshold
    }

    statistics_table = pd.DataFrame(list(statistics.items()),
                                    columns = ['Exclusion Reason', 'Sample Count'])
    
    statistics_table["Percentage"] = statistics_table["Sample Count"].apply(lambda x: '{:.1f} %'.format(100 * x / n_qc_pass_samples))

    st.write(f"{included_samples} samples ({included_samples / n_qc_pass_samples * 100:.1f} %) are available for analysis out of {n_qc_pass_samples} QC pass samples for this gene, after {excluded_samples} samples ({excluded_samples / n_qc_pass_samples * 100:.1f} %) have been excluded for the reasons listed below. ")
    st.dataframe(statistics_table, use_container_width = True, hide_index = True)

    return

def _config_download_data_section(gene_id_selected, df_haplotypes, df_join):
    st.markdown("Hover over the download buttons for more information on the data.")
    
    st.download_button(
        "Download population-level summary",
        df_haplotypes.reset_index(drop = True).to_csv().encode('utf-8'),
        file_name = f'pf-haploatlas-{gene_id_selected}_population_summary.csv',
        help = '''Explanation of columns: "number_of_mutations" describes number of mutations relative to 3D7; "ns_changes" describes the amino acid changes of each unique haplotype; "SA", "AF-W", "AF-C", etc. shows number of samples observed with that haplotype in each geographic distribution (see sidebar for details); "Total" is the total number of samples with that haplotype; "ns_changes_list" is a list of amino acid changes of the haplotype; "sample_names" describes which lab strains the haplotype is found in''',
        use_container_width = True)
    
    st.download_button(
        "Download sample-level summary",
        df_join.to_csv().encode('utf-8'),
        file_name = f'pf-haploatlas-{gene_id_selected}_sample_summary.csv',
        help = '''Explanation of columns: "Sample" is the sample name; "Study" is the clinical study of origin; "Country" of sample collection; "Admin level"	is the location of sample collection; "latitude", "longitude, "Year" of sample collection; "ENA" is the ID in the European Nucleotide Archive; "All samples same case" is reformatted sample name, "Population" refers to geographic distribution (see sidebar for details); "% callable" of SNPs" describes what percentage of QC-passed SNP positions were callable; "QC pass" is whether the sample passed quality control for Pf8; "Exclusion reason" describes the reason for a sample's removal from the Pf8 QC-pass cohort or whether it was part of the QC-passed analysis set ("Analysis_set"); "Sample type" describes the type of sample put through to sequencing; "Sample was in Pf7" is whether the sample was in the previous Pf7 data resource; "ns_changes" describes the amino acid changes of the sample for the gene selected; "HaploAtlas exclusion reason" describes for what reason the sample was excluded from analysis on the app for the sample number threshold used at the time of downloading the data ("Missing_genotype" means a genotype call was missing so amino acid haplotype could not be reliably identified; "Het_calls" means there were SNPs which were called as heterozygous so amino acid haplotype could not be reliably identified; "Unverified_identity" means there was missing metadata in Pf8; "Stop_codon" implies SNPs suggested that the protein would be truncated by a stop codon so amino acid haplotype could not be reliably identified)''',
        use_container_width = True)
    
    return
    
def _config_plot_settings_section():
    sample_count_mode = st.radio("Select y-axis mode for Haplotype UpSet plot:", 
                                 ["Sample counts", "Sample counts on a log scale"],
                                 index = 0)
    return sample_count_mode