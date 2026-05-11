import streamlit as st
import json, os, lzma, pickle, collections, io
import pandas as pd

DATAPACK_BASE_PATH = "app/files/pf9-datapack-2026-02-19"

@st.cache_data
def _cache_load_utility_mappers(base_path = DATAPACK_BASE_PATH):
    """
    Loads various useful dictionaries and lists related to handling gene IDs and converting
    back and forth between gene IDs and gene names etc. 
    Caches the objects when first loaded
    """
    
    with open(f"{base_path}/auxiliary/core_genes.json", "r") as f:
            gene_mapper = json.load(f)
    
    files_to_gene_ids = {f: f.split(".")[0] for f in os.listdir(f"{base_path}/genes")}
    
    gene_ids_to_files = dict( zip(files_to_gene_ids.values(), files_to_gene_ids.keys()) )
    
    gene_ids_to_gene_names = {
        gene: (f'{gene} - {gene_mapper[gene]}' if gene_mapper[gene] != "" else gene)
            for gene in sorted(files_to_gene_ids.values())
    }
    
    gene_names_to_gene_ids = dict(zip(gene_ids_to_gene_names.values(), gene_ids_to_gene_names.keys()))
    
    gene_ids = [gene_id for gene_id, gene_name in gene_ids_to_gene_names.items()] # before core genes identified: if _is_core_genome(gene_name)
    
    return {
        "gene_ids_to_files"     : gene_ids_to_files,
        "gene_ids_to_gene_names": gene_ids_to_gene_names,
        "gene_names_to_gene_ids": gene_names_to_gene_ids,
        "gene_ids"              : gene_ids
    }

@st.cache_data
def _cache_load_sample_metadata():
    return pd.read_csv(f'{DATAPACK_BASE_PATH}/auxiliary/2026_05_11_Pf9_sample_metadata.tsv.gz', sep = "\t", low_memory = False)

@st.cache_data
def cache_load_gene_summary(filename: str, base_path = DATAPACK_BASE_PATH):
    """Loads the relevant gene summary file based on provided file path. Caches the objects when first loaded"""    

    with lzma.open(f'{base_path}/genes/{filename}', 'rb') as file:
        df_haplotypes, df_join, background_ns_changes = pickle.load(file)
    
    df_join = df_join.rename(columns = {"Exclusion reason": "HaploAtlas exclusion reason"})
    df_join = pd.concat([_cache_load_sample_metadata(), df_join], axis = 1)

    df_join["Exclusion reason"] = df_join["Exclusion reason"].fillna("Analysis_set")
    df_join["HaploAtlas exclusion reason"] = df_join["HaploAtlas exclusion reason"].fillna("Analysis_set")
    
    return df_haplotypes, df_join, background_ns_changes

@st.cache_data
def cache_load_population_colours():
    """Pf9 population colour palette. Caches the objects when first loaded"""
    population_colours = collections.OrderedDict()
    population_colours["LA-W"]    = "#b8e186"
    population_colours["LA-E"]    = "#4dac26"
    population_colours["AF-W"]    = "#e31a1c"
    population_colours["AF-C"]    = "#fd8d3c"
    population_colours["AF-NE"]   = "#bb8129"
    population_colours["AF-E"]    = "#fecc5c"
    population_colours["AS-S-E"]  = "#dfc0eb"
    population_colours["AS-S-FE"] = "#984ea3"
    population_colours["AS-SE-W"] = "#9ecae1"
    population_colours["AS-SE-E"] = "#3182bd"
    population_colours["AS-SE-M"] = "#02818a"
    population_colours["OC-NG"]   = "#f781bf"
    population_colours["EU"]      = "#003399"
    
    return population_colours

def _cache_load_changelog():
    with open("app/files/changelog.md", "r") as f:
        return f.read()

def generate_download_buttons(fig, gene_id_selected, height, width, plot_number):
    """Generates download buttons for different image formats (PDF, PNG, SVG) for a given plot."""

    plot_name_dictionary = {
        1: "haplotype_upset_plot",
        2: "abacus_plot",
        3: "worldmap_plot"
    }

    plot_name = plot_name_dictionary[plot_number]

    # Create in-memory buffers to download the plot
    buffers = {}
    formats = ["pdf", "png", "svg"]
    for format in formats:
        buffer = io.BytesIO()
        fig.write_image(file=buffer, format=format, height=height, width=width)
        buffers[format] = buffer

    figure_name = f"{gene_id_selected}_{plot_name}"

    col, *button_cols = st.columns([7, 1, 1, 1])
    col.markdown("<p style='text-align: right; line-height: 40px;'>Download the figure:</p>", unsafe_allow_html=True)
    
    for col, format in zip(button_cols, formats):
        col.download_button(
            label     = format.upper(),
            data      = buffers[format],
            file_name = f"pf-haploatlas-{figure_name}.{format}",
            mime      = f"image/{format}",
            type      = "primary",
            key       = f'{plot_number}_{format}'
        )

def haplotype_selection_toast(ns_changes):

    toast_message = f"Haplotype {ns_changes} selected! Scroll down to see the Abacus and world map plots"

    if "ns_changes" not in st.session_state:
        st.session_state["ns_changes"] = ns_changes
        st.toast(toast_message, icon = "⬇️")

    elif st.session_state["ns_changes"] != ns_changes:
        st.toast(toast_message, icon = "⬇️")
        st.session_state["ns_changes"] = ns_changes

def _st_justify_markdown_html(text: str, location = None):

    justify_css = """
    <style>
        .justify-text {
            text-align: justify;
        }
    </style>
    """


    text = f"""
<div class="justify-text">
{text}
</div>
"""

    if location is None:
        st.markdown(
            justify_css + text, 
            unsafe_allow_html = True
        )
    
    else:
        location.markdown(
            justify_css + text, 
            unsafe_allow_html = True
        )
    return

@st.dialog("Pf9 is here! 🥳")
def _present_cookie_banner():
    tab1, tab2 = st.tabs(["Welcome", "Cookie details"])

    tab1.markdown("We now use [**MalariaGEN Pf8**](https://wellcomeopenresearch.org/articles/10-325), which means that the Pf-HaploAtlas now contains data from **40,074** QC-passed samples. That's _15,665_ more, or over _64%_ more than the last version. Happy haplotype hunting!")

    tab1.markdown("""<p style="text-align: right;"><i>- the HaploAtlas team</i></p>""", unsafe_allow_html = True)

    tab1.divider()
    
    tab1.write("We'd also like to use analytics cookies to help improve our page.")

    _, tab1col1, tab1col2, _  = tab1.columns([1, 2, 2, 1])
    if tab1col1.button("Accept", use_container_width = True, key = "tab1col1"):
        st.session_state["cookies_accepted"] = True
        st.rerun()
    
    if tab1col2.button("Reject", use_container_width = True, key = "tab1col2"):
        st.session_state["cookies_accepted"] = False
        st.rerun()

    with tab2:
        _st_justify_markdown_html("""
## Cookie Statement:

This site uses cookies to store information on your computer.

### What are cookies?
Cookies are small data files that are sent from a website you visit to your browser and stored on the device you use to access the internet.

### Why do we use cookies?
Cookies perform a variety of functions that enable our websites to operate efficiently, help us improve their performance, and provide the best user experience.

### What cookies do we use?

Session Cookies: These cookies allow us to understand how you've used our website during a browsing session. Session cookies expire once your browsing session ends and are not stored beyond this point.

Google Analytics Cookies: We use Google Analytics to help us understand how visitors navigate our websites. Google Analytics uses analytics cookies to collect information and report website usage statistics without personally identifying individual visitors. The information collected is stored for 26 months and includes:
- Number of users to our website
- The country you are accessing the website from
- Details of your browser and device
- User behavior
- The website which directed you to our website

Learn more about who we are, how you can contact us, and how we process personal data in the [MalariaGEN website's privacy policy](https://www.malariagen.net/privacy-policy/).
""")
        
        _, tab2col1, tab2col2, _  = tab2.columns([1, 2, 2, 1])
        if tab2col1.button("Accept", use_container_width = True, key = "tab2col1"):
            st.session_state["cookies_accepted"] = True
            st.rerun()

        if tab2col2.button("Reject", use_container_width = True, key = "tab2col2"):
            st.rerun()

def _show_cookie_banner_upon_visit():
    if "banner_shown" not in st.session_state:
        _present_cookie_banner()
        st.session_state["banner_shown"] = True
    else:
        pass

def present_changelog():
    with st.expander("Click to see change log"):
        logs = _cache_load_changelog()
        st.write(logs)
