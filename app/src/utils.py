import streamlit as st
import json, os, lzma, pickle

base_path = "app/files/06-12-23_smaller_pkls"

def _is_core_genome(filename: str):
        if "VAR" in filename:
            return False
        elif "SURF" in filename:
            return False
        elif "RIF" in filename:
            return False
        else:
            return True

@st.cache_data
def _cache_load_utility_mappers(base_path = base_path):
    
    with open("app/files/gene_mapping.json", "r") as f:
            gene_mapper = json.load(f)
    
    files_to_gene_ids = {f: f.split(".")[0] for f in os.listdir(base_path) if f.endswith("pkl.xz")}
    
    gene_ids_to_files = dict( zip(files_to_gene_ids.values(), files_to_gene_ids.keys()) )
    
    gene_ids_to_gene_names = {
        gene: (f'{gene} - {gene_mapper[gene]}' if not gene_mapper[gene] is "." else gene)
            for gene in sorted(files_to_gene_ids.values())
    }
    
    gene_names_to_gene_ids = dict(zip(gene_ids_to_gene_names.values(), gene_ids_to_gene_names.keys()))
    
    gene_ids = [gene_id for gene_id, gene_name in gene_ids_to_gene_names.items() if _is_core_genome(gene_name)]
    
    return {
        "gene_ids_to_files": gene_ids_to_files,
        "gene_ids_to_gene_names": gene_ids_to_gene_names,
        "gene_names_to_gene_ids": gene_names_to_gene_ids,
        "gene_ids": gene_ids
    }

@st.cache_data
def cache_load_gene_summary(filename: str, base_path = base_path):
    with lzma.open(f'{base_path}/{filename}', 'rb') as file:
        loaded_plot_data = pickle.load(file)
    return loaded_plot_data
