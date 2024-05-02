import streamlit as st
from PIL import Image
import json, os

from src.utils import _cache_load_utility_mappers, _cache_load_pf7_metadata

def set_up_interface():
    """Main function called in main.py to set up basic page settings, introduction and sidebar"""
    
    st.set_page_config(
        page_title = "Pf7 mutation discovery app",
        page_icon = Image.open("app/files/logo.png"),
        initial_sidebar_state = "expanded"
    )
    
    st.title('Pf7.0 Haplotype Explorer (PfHEx7.0)')
    st.subheader("Haplotype summaries for *Plasmodium falciparum* genes across time and space")

    _cache_load_pf7_metadata() # running it here to prevent it from running when new gene selected
    
    st.divider()
    
    placeholder = st.empty()
    placeholder.markdown("### Search for a gene below to get started.")
    
    with st.sidebar:
        st.title("**Further Information**")
        
        st.divider()
        
        st.header("**Samples**")
        st.markdown("""

The Mutation Discovery App uses 16,203 QC pass samples from the [Pf7 dataset.](https://wellcomeopenresearch.org/articles/8-22/v1)

""")
        
        st.header("**Genes**")
        st.markdown("""

The Mutation Discovery App uses 5102 genes located within the core regions of 3D7 v3 reference genome (available at this [ftp server](ftp://ngs.sanger.ac.uk/production/malaria/Resource/34/Pfalciparum.genome.fasta)). All genes have a unique identifier, e.g. **PF3D7_1343700**, and in some cases a gene name, e.g. **MDR1**.

""")
        st.header("**Subpopulations**")
        st.markdown("""
Countries are grouped into ten major sub-populations based on their geographic and genetic characteristics as defined in the [Pf7 paper](https://wellcomeopenresearch.org/articles/8-22/v1). These are colour-coded for easy interpretation. 
                    """)
        st.header("**Plots**")
        st.markdown("""

The app generates five plots per gene:

**1. Sample counts per haplotype** - cumulative counts of samples per haplotype. Toggle the y-axis between raw values or a log scale.

**2. Population proportions per haplotype** - relative contributions of each of the major Pf7 sub-populations to each haplotype.

**3. UpSet plot** - view the mutation makeup of each haplotype.

**4. Abacus plot** - click on a haplotype to view its frequency across first-level administrative divisions and years with at least 25 samples.

The shade of the point represents the haplotype frequency from white (0%) to black (100%). Where frequency is exactly 0% or 100% the point is marked with a cross to represent fixation.

**5. World map plot** - use the slider to specify a range of years for visualizing haplotype frequencies across countries on a world map with at least 25 samples. For the explanation of points, please refer to Abacus plot.                  
""")
    
    return placeholder

def file_selector(placeholder):
    """Main function called in main.py to allow for user's gene selection and handle the app's URL"""
    utility_mappers = _cache_load_utility_mappers()
    # Check if 'gene_id' key exists in query parameters
    if 'gene_id' in st.query_params:
        gene_id_extracted = st.query_params['gene_id']
    else:
        gene_id_extracted = "--"
    
    gene_id_extracted = utility_mappers["gene_ids_to_gene_names"].get(gene_id_extracted, "--")

    if gene_id_extracted and "gene_id" not in st.session_state:
        st.session_state["gene_id"] = gene_id_extracted
    
    gene_id_selected = st.selectbox("",
                                    ["--"] + [utility_mappers["gene_ids_to_gene_names"][gene_id]
                                              for gene_id in utility_mappers["gene_ids"] 
                                              if gene_id in utility_mappers["gene_ids_to_gene_names"].keys()],
                                    key="gene_id"
                                   )
    
    if gene_id_selected == "--":
        placeholder.markdown("### Search for a gene below to get started.")
        st.query_params.get_all('gene_id')
        st.stop()
    
    gene_id_selected = utility_mappers["gene_names_to_gene_ids"].get(gene_id_selected, "--")
    
    if gene_id_selected != gene_id_extracted:
        st.query_params["gene_id"] = gene_id_selected

    filename = utility_mappers["gene_ids_to_files"].get(gene_id_selected, None)
    
    if filename is None:
        st.warning(f"No file found for gene ID: {gene_id_selected}")
        st.stop()

    placeholder.empty()
    return filename, gene_id_selected