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

def process_configs_menu_kelch():
    """Main function called in main.py to handle user config settings in the expander"""
    col1, col2 = st.columns([7, 5])
    
    # letter-codes are used for the datasets. 
    # we have 13 combinations, we will load a file for each combination. 
    dataset_files = {
    "Pf7": "a",
    "Amplicon in-house": "b",
    "Amplicon in-country": "c",
    "Literature search": "d",
    }

    with col1.expander("Click to edit plot settings"):
        min_samples = st.number_input("Minimum number of samples per haplotype for analysis", min_value = 1, value = 25)
        sample_count_mode = st.radio("Select a mode for Plot 1:", ["Sample counts", "Sample counts on a log scale"], index=0)
        select_dataset = st.multiselect("Select the dataset(s)", list(dataset_files.keys()), default=list(dataset_files.keys()))
        # joins the letters to make a alphebetically ordered word
        select_dataset_value = "".join(sorted(dataset_files[dataset] for dataset in select_dataset))
        #print(select_dataset_value)
        #select_dataset_value= dataset_files[selected_word]
        selected_gene_plasmodb_url = f'https://plasmodb.org/plasmo/app/record/gene/PF3D7_1343800'
        
    
    col2.markdown(
        f'''<a href="{selected_gene_plasmodb_url}" style="display: inline-block;
            padding: 11px 10px; background-color: #b00023;
            color: white;
            text-align: center;
            text-decoration: none;
            font-size: 16px; border-radius:
        4px;">Browse Gene on PlasmoDB</a>''',
        unsafe_allow_html=True
    )
    
    return min_samples, sample_count_mode, select_dataset_value