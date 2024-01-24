import streamlit as st

def process_configs_menu(gene_id_selected):
    """Main function called in main.py to handle user config settings in the expander"""
    col1, col2 = st.columns([3, 2])
    
    with col1.expander("Click to edit plot settings"):
        min_samples = st.number_input("Minimum number of samples per haplotype for analysis", min_value = 1, value = 25)
        sample_count_mode = st.radio("Select a mode for Plot 1:", ["Sample counts", "Sample counts on a log scale"], index=0)

    selected_gene_plasmodb_url = f'https://plasmodb.org/plasmo/app/record/gene/{gene_id_selected}'
    
    col2.markdown(
        f'''<a href="{selected_gene_plasmodb_url}" style="display: inline-block;
            padding: 11px 20px; background-color: #b00023;
            color: white;
            text-align: center;
            text-decoration: none;
            font-size: 16px; border-radius:
        4px;">Browse Gene on PlasmoDB</a>''',
        unsafe_allow_html=True
    )
    
    return min_samples, sample_count_mode
