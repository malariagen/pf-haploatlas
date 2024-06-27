import streamlit as st
from base64 import b64encode

from src.utils import _cache_load_utility_mappers, _cache_load_pf7_metadata, _st_justify_markdown_html



def set_up_interface():
    """Main function called in main.py to set up basic page settings, introduction and sidebar"""

    st.set_page_config(
        page_title = "Pf-HaploAtlas",
        layout = "centered",
        initial_sidebar_state = "expanded"
    )
    
    st.markdown(
        """
        <div style="text-align:center; font-size:4em; font-weight:bold;">
            Pf-HaploAtlas
        </div>
        """, 
        unsafe_allow_html=True
    )

    _cache_load_pf7_metadata() # running it here to prevent it from running when new gene selected
    
    st.divider()
    
    placeholder = st.empty()

    _st_justify_markdown_html("""
### Introduction

The _Plasmodium falciparum_ Haplotype Atlas (or Pf-HaploAtlas) allows anyone with an internet connection to study and track genetic mutations across any gene in the _P. falciparum_ genome! The app provides visualisations of haplotypes for all 5,102 core genes by using data from 16,203 samples, from 33 countries, and spread between the years 1984 and 2018, facilitating comprehensive spatial and temporal analyses of genes and variants of interest. Please check out our tutorial video in the sidebar to learn how to use the app. This web app was primarily developed for use on a desktop browser. If you would like support for mobile, please request this feature in the feedback form in the sidebar! We also encourage users to access and share the app using the following stable link to prevent outages in service: https://apps.malariagen.net/pf-haploatlas.

The Pf-HaploAtlas journal manuscript can be found in the sidebar. Pf-HaploAtlas currently uses data generated using the [MalariaGEN Pf7 whole genome sequencing data release](https://wellcomeopenresearch.org/articles/8-22/v1), and will expand with each new MalariaGEN _Plasmodium_ data release. 

#### Search for a gene below to get started.

If you're new here, try clicking below and typing "AAT1"! Alternatively, choose from the key drug resistance genes we've placed at the top of the list (DHFR-TS, MDR1, CRT, PPPK-DHPS, Kelch13).

""", location = placeholder)

    _set_up_sidebar()
    
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
    
    priority_gene_ids = [
        "PF3D7_0417200", "PF3D7_0523000", "PF3D7_0709000", "PF3D7_0810800", "PF3D7_1343700"
    ]
    priority_gene_names = [utility_mappers["gene_ids_to_gene_names"][gene_id] for gene_id in priority_gene_ids]
    
    gene_id_selected = st.selectbox(" ",
                                    ["--"] + priority_gene_names + ["--"] + [utility_mappers["gene_ids_to_gene_names"][gene_id]
                                              for gene_id in utility_mappers["gene_ids"] 
                                              if gene_id in utility_mappers["gene_ids_to_gene_names"].keys()],
                                    key="gene_id", label_visibility='collapsed'
                                   )
    
    if "--" in gene_id_selected:
        # placeholder.markdown("### Search for a gene below to get started.")
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

def _show_images_with_urls(filepaths, urls, widths, heights):
    images_html = "<div style='display: flex; justify-content: center; align-items: flex-end; text-align: center;'>"
    for filepath, url, width, height in zip(filepaths, urls, widths, heights):
        image_data = b64encode(open(filepath, "rb").read()).decode()
        images_html += f"""
            <div style="margin: 10px;">
                <a href="{url}">
                    <img src="data:image/png;base64,{image_data}" style="width: {width}%; height: {height}%; object-fit: contain;">
                </a>
            </div>"""
    images_html += "</div>"
    st.markdown(images_html, unsafe_allow_html=True)

def _set_up_sidebar():
    # This fixes the width of the sidebar to a specified number of pixels
    st.markdown("""
<style>
    section[data-testid="stSidebar"] {
        width: 350px !important;
    }
</style>
    """, unsafe_allow_html=True)

    with st.sidebar:
                
        st.title("How to use")

        st.video("https://www.youtube.com/watch?v=XLGL1xqZEQk")
        st.divider()

        _st_justify_markdown_html("""
## Overview of plots:

**1. Haplotype UpSet plot** - for each haplotype of your chosen gene, shows the number of samples with that haplotype, the geographic distribution of these samples, and the mutation make-up.

Clicking on a haplotype will generate the two following plots: 

**2. Abacus plot** - for each location, shows the proportion of samples containing your chosen haplotype each year

**3. World map** - for each country, shows the proportion of samples with your chosen haplotype over your selected time period
""")

        st.divider()

        _st_justify_markdown_html("""
## Geographic distribution

The locations of where samples were collected are grouped into ten major "sub-populations" based on their geographic and genetic characteristics, as defined in the <a href="https://wellcomeopenresearch.org/articles/8-22/v1" target="_blank">Pf7 paper</a>. These are colour-coded as follows:
""")

        st.markdown("""
<ul style="list-style-type:none;">
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#f781bf; border-radius:50%;"></span> OC-NG - Oceania, New Guinea</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#3182bd; border-radius:50%;"></span> AS-SE-E - South-East Asia (East)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#9ecae1; border-radius:50%;"></span> AS-SE-W - South-East Asia (West)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#984ea3; border-radius:50%;"></span> AS-S-FE - South Asia (Far East)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#dfc0eb; border-radius:50%;"></span> AS-S-E - South Asia (East)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#fecc5c; border-radius:50%;"></span> AF-E - Africa (East)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#bb8129; border-radius:50%;"></span> AF-NE - Africa (North-East)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#fd8d3c; border-radius:50%;"></span> AF-C - Africa (Central)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#e31a1c; border-radius:50%;"></span> AF-W - Africa (West)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#4daf4a; border-radius:50%;"></span> SA - South America</li>
</ul>
""", unsafe_allow_html=True)

        _st_justify_markdown_html("""
On the x-axis of the geographic distribution subplot of the Haplotype UpSet plot, you will also see the names of lab strains which are also of that haplotype (e.g., 3D7, 7G8, Dd2, IT, GB4, HB3). 
""")
        st.divider()

        _st_justify_markdown_html("""
## How to cite

When publishing work that uses data and/or plots from the Pf-HaploAtlas, please cite the following:

```
Lee C, Unlu E, White NFD et al. Pf-HaploAtlas: an interactive web app for spatiotemporal analysis of P. falciparum genes. BioRxiv 123456 [Preprint]. January 01, 2024 [cited 2030 Jan 01]. Available from: https://doi.org/10.1234/123456.
```

""")

        st.divider()

        _st_justify_markdown_html("""
## Acknowledgements

Pf-HaploAtlas currently uses data generated using the [Pf7 data release](https://wellcomeopenresearch.org/articles/8-22/v1) from the MalariaGEN _Plasmodium falciparum_ Community Project, which was made possible by clinical parasite samples contributed by partner studies, whose investigators are represented in the data release's author list.

""")

        st.divider()

        _st_justify_markdown_html("""
## Contact us
If you'd like to report a bug, request a feature, or give us feedback, check out the following!

- [our GitHub page](https://github.com/malariagen/pf-haploatlas/issues)
- [this Google Forms](https://forms.gle/mDwYr2cPL37dDzPs6)
""")

        st.divider()

        _show_images_with_urls(
            ["app/files/logo_malariagen.png"],
            ["https://www.malariagen.net/"],
            [50],
            [100]
        )

        _show_images_with_urls(
            ["app/files/logo_gsu.png", "app/files/logo_sanger.png"],
            ["https://www.sanger.ac.uk/collaboration/genomic-surveillance-unit/", "https://www.sanger.ac.uk/"],
            [110, 70],
            [110, 70]
        )

        st.divider()

        st.markdown(
            """
<div style="text-align:center">
    Copyright © 2021 - 2024 Genome Research Ltd.
</div>
            """
        , unsafe_allow_html = True)
