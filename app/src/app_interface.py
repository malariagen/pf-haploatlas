import streamlit as st
import json, os

from PIL import Image

from src.utils import _cache_load_utility_mappers

def set_up_interface():
    """Main function called in main.py to set up basic page settings, introduction and sidebar"""
    
    st.set_page_config(
        page_title = "Pf7 mutation discovery app",
        page_icon = Image.open("app/files/logo.png"),
        initial_sidebar_state = "expanded"
    )
    
    st.title('Pf7.0 Haplotype Explorer (PfHEx7.0)')
    st.subheader("Haplotype summaries for *Plasmodium falciparum* genes across time and space")
    
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

The Mutation Discovery App uses **XXX** genes, according to the **XXX** version of the 3D7 Pf genome from *PlasmoDB link?* All genes have a unique identifier, e.g. **PF3D7_XXXXXXX**, and in some cases a gene name, e.g. **MDR1**.

""")
        
        st.header("**Plots**")
        st.markdown("""

The app generates four plots per gene:

**1. Sample counts per haplotype** - cumulative counts of samples per haplotype. Toggle the y-axis between raw values or a log scale.

**2. Population proportions per haplotype** - relative contributions of each of the major Pf7 sub-populations to each haplotype.

**3. UpSet plot** - view the mutation makeup of each haplotype.

**4. Abacus plot** - click on a haplotype to view its frequency across first-level administrative divisions and years with at least 25 samples.

The shade of the point represents the haplotype frequency from white (0%) to black (100%). Where frequency is exactly 0% or 100% the point is marked with a cross to represent fixation.

""")
    
    return placeholder

def file_selector(placeholder):
    """Main function called in main.py to allow for user's gene selection and handle the app's URL"""
    utility_mappers = _cache_load_utility_mappers()
    gene_id_extracted = st.query_params.get_all('gene_id')

    if placeholder not in utility_mappers["gene_ids_to_gene_names"]:
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
        #st.experimental_set_query_params(**current_query_params)

    filename = utility_mappers["gene_ids_to_files"].get(gene_id_selected, None)
    
    if filename is None:
        st.warning(f"No file found for gene ID: {gene_id_selected}")
        st.stop()

    placeholder.empty()
    return filename, gene_id_selected

def set_up_interface_kelch():
    """Main function called in main.py to set up basic page settings, introduction and sidebar"""
    
    st.set_page_config(
        page_title = "Kelch13 Haplotype Explorer",
        page_icon = Image.open("app/files/logo.png"),
        initial_sidebar_state = "expanded"
    )
    
    st.title('Kelch13 Haplotype Explorer')
    st.subheader("Haplotype summary for *Plasmodium falciparum* K13 gene across time and space")
    
    st.divider()
   
       
    placeholder = st.empty()
    placeholder.markdown("### Search for a gene below to get started.")

    with st.sidebar:
        st.title("**Further Information**")
        
        st.divider()
        
        st.header("**Samples**")
        st.markdown("""

The Kelch13 Haplotype Explorer contains data from 108,221 *Plasmodium falciparum* samples

""")
        
        st.header("**Data Sources**")
        st.markdown("""

Data is collected from one of three sources:

1. Amplicon sequencing of K13 from the GSU (n = 23,138)

2. Whole genome sequencing data from the Pf7 release (n = 9,593)

3. A literature review. Contains data from the [WWARN database](http://www.wwarn.org/molecular/surveyor/k13/index.html?t=201608031200#0) covering years 1980 - 2021, and a supplementary search from years 2020 - 2024 (n = 75,490).

You may view data from individual or combined sources in the app.

""")
        
        
        st.header("**Plots**")
        st.markdown("""

The app generates four plots per gene:

**1. Sample counts per haplotype** - cumulative counts of samples per haplotype. Toggle the y-axis between raw values or a log scale.

**2. Population proportions per haplotype** - relative contributions of each of the major Pf7 sub-populations to each haplotype.

**3. UpSet plot** - view the mutation makeup of each haplotype. Status of K13 mutations relating to artemisinin resistance are colour coded according to WHO classifications. 

**4. Abacus plot** - click on a haplotype to view its frequency across countries and years.

The shade of the point represents the haplotype frequency from white (0%) to black (100%). Where frequency is exactly 0% or 100% the point is marked with a cross to represent fixation.

""")
    
    return placeholder