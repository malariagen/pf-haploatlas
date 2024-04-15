import streamlit as st
import json, os, lzma, pickle, collections
import pandas as pd

base_path = "app/files/2024-03-13_pkl_files"

@st.cache_data
def _cache_load_utility_mappers(base_path = base_path):
    """
    Loads various useful dictionaries and lists related to handling gene IDs and converting
    back and forth between gene IDs and gene names etc. 
    Caches the objects when first loaded
    """
    
    with open("app/files/core_genes.json", "r") as f:
            gene_mapper = json.load(f)
    
    files_to_gene_ids = {f: f.split(".")[0] for f in os.listdir(base_path) if f.endswith("pkl.xz")}
    
    gene_ids_to_files = dict( zip(files_to_gene_ids.values(), files_to_gene_ids.keys()) )
    
    gene_ids_to_gene_names = {
        gene: (f'{gene} - {gene_mapper[gene]}' if not gene_mapper[gene] is "." else gene)
            for gene in sorted(files_to_gene_ids.values())
    }
    
    gene_names_to_gene_ids = dict(zip(gene_ids_to_gene_names.values(), gene_ids_to_gene_names.keys()))
    
    gene_ids = [gene_id for gene_id, gene_name in gene_ids_to_gene_names.items()] # before core genes identified: if _is_core_genome(gene_name)
    
    return {
        "gene_ids_to_files": gene_ids_to_files,
        "gene_ids_to_gene_names": gene_ids_to_gene_names,
        "gene_names_to_gene_ids": gene_names_to_gene_ids,
        "gene_ids": gene_ids
    }

@st.cache_data
def _cache_load_pf7_metadata():
    pf7_metadata = pd.read_excel('app/files/Pf7_metadata.xlsx')
    return pf7_metadata

@st.cache_data
def cache_load_gene_summary(filename: str, base_path = base_path):
    """Loads the relevant gene summary file based on provided file path. Caches the objects when first loaded"""    
    with lzma.open(f'{base_path}/{filename}', 'rb') as file:
        loaded_plot_data = pickle.load(file)
    df_haplotypes, df_join, background_ns_changes, _ = loaded_plot_data
    df_join = pd.concat([df_join.reset_index(), _cache_load_pf7_metadata()], axis=1)
    return df_haplotypes, df_join, background_ns_changes

@st.cache_data
def cache_load_population_colours():
    """Pf7 population colour palette. Caches the objects when first loaded"""
    population_colours = collections.OrderedDict()
    population_colours['SA'] = "#4daf4a"
    population_colours['AF-W']= "#e31a1c"
    population_colours['AF-C'] = "#fd8d3c" 
    population_colours['AF-NE'] = "#bb8129" 
    population_colours['AF-E'] = "#fecc5c"
    population_colours['AS-S-E'] = "#dfc0eb" 
    population_colours['AS-S-FE'] = "#984ea3" 
    population_colours['AS-SE-W'] = "#9ecae1"
    population_colours['AS-SE-E'] = "#3182bd"
    population_colours['OC-NG'] = "#f781bf"
    
    return population_colours

@st.cache_data
def cache_load_population_colours_kelch():
    """Pf7 population colour palette. Caches the objects when first loaded"""
    population_colours = collections.OrderedDict()
    population_colours['SA'] = "#4daf4a"
    population_colours['AF-W']= "#e31a1c"
    population_colours['AF-C'] = "#fd8d3c"
    population_colours['AF-N'] = "#654321"
    population_colours['AF-NE'] = "#bb8129" 
    population_colours['AF-E'] = "#fecc5c"
    population_colours['AF-S'] = "#800000"
    population_colours['AS-W'] = "#D3D3D3"
    population_colours['AS-S-E'] = "#dfc0eb" 
    population_colours['AS-S-FE'] = "#984ea3" 
    population_colours['AS-SE-W'] = "#9ecae1"
    population_colours['AS-SE-E'] = "#3182bd"
    population_colours['OC-NG'] = "#f781bf"
    
    return population_colours

@st.cache_data
def cache_load_gene_summary_kelch(select_dataset_value, base_path = base_path):
    """Loads the relevant gene summary file based on provided file path. Caches the objects when first loaded"""
    with lzma.open(f'app/files/Kelch_all_data.pkl.xz', 'rb') as file:
        loaded_plot_data = pickle.load(file)
    return loaded_plot_data[select_dataset_value]