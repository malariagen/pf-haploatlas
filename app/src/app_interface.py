import streamlit as st
from base64 import b64encode

from src.utils import _cache_load_utility_mappers, _cache_load_sample_metadata, _st_justify_markdown_html, _show_cookie_banner_upon_visit, present_changelog
from streamlit_gtag import st_gtag

def set_up_interface():
    """Main function called in main.py to set up basic page settings, introduction and sidebar"""

    st.set_page_config(
        page_title            = "Pv-HaploAtlas",
        layout                = "centered",
        page_icon             = "app/files/favicon.svg",
        initial_sidebar_state = "collapsed",
    )

    if "cookies_accepted" in st.session_state:
        if "gtag_injected" not in st.session_state:
            st_gtag(
                key="gtag_send_event_a",
                id="G-4XZZ9XXZ21",
                event_name="cookies_accepted",
                params={
                    "event_category": "test_category_a",
                    "event_label": "test_label_a",
                    "value": "test",
                },
            )
            st.session_state["gtag_injected"] = True

    hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
    st.markdown(hide_streamlit_style, unsafe_allow_html=True) 
    
    st.markdown(
        """
        <div style="text-align:center; font-size:4em; font-weight:bold;">
            Pv-HaploAtlas
        </div>
        """, 
        unsafe_allow_html=True
    )

    _cache_load_sample_metadata() # running it here to prevent it from running when new gene selected
    
    st.divider()
    
    placeholder = st.empty()

    _st_justify_markdown_html("""
### Introduction

The _Plasmodium vivax_ Haplotype Atlas (or Pv-HaploAtlas) allows anyone with an internet connection to study and track genetic mutations across any gene in the _P. vivax_ genome! The app provides visualisations of haplotypes for all 4,936 core genes by using data from 2,298 samples, from 29 countries, and spread between the years 1971 and 2019, facilitating comprehensive spatial and temporal analyses of genes and variants of interest. 
                              
The Pv-HaploAtlas has been designed to be user-friendly and intuitive, so that you can <b><u>learn how to use the app by simply continuing to read this page</b></u>. If you prefer, refer to our tutorial video or manuscript in the sidebar for a quick demo. To request a feature, get involved on our GitHub Issues or feedback form in the sidebar. 

Pv-HaploAtlas currently uses data generated using the [MalariaGEN Pv5 whole genome sequencing data release](insert link), and will expand with each new MalariaGEN _Plasmodium vivax_ data release. 

#### Search for a gene below to get started.

If you're new here, try clicking below and typing "CRT"! Alternatively, choose from the key drug resistance genes we've placed at the top of the list (CRT, DHFR, DHPS, MDR1).

""", location = placeholder)

    _set_up_sidebar()

    _show_cookie_banner_upon_visit()
    
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
        "PVP01_0109300", "PVP01_0526600", "PVP01_1429500", "PVP01_1010900"
    ]
    priority_gene_names = [utility_mappers["gene_ids_to_gene_names"][gene_id] for gene_id in priority_gene_ids]
    
    gene_id_selected = st.selectbox(" ",
                                    ["--"] + priority_gene_names + ["--"] + [utility_mappers["gene_ids_to_gene_names"][gene_id]
                                              for gene_id in utility_mappers["gene_ids"] 
                                              if gene_id in utility_mappers["gene_ids_to_gene_names"].keys()],
                                    key = "gene_id",
                                    help = """This list of core genes was created using "protein coding genes" (as defined by the GFF of PlasmoDB version 68 for _P. vivax_ reference strain P01) based on core genome region annotations from [REF]. Gene IDs are accompanied by values from the GFF's "ID" field value or if unavailable, from the "description" field (in which case it is enclosed by quotation marks). """
                                   )
    
    if "--" in gene_id_selected:
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

        st.video("https://youtu.be/48f4r2frcdk")

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

The locations of where samples were collected are grouped into 7 major "sub-populations" based on their geographic and genetic characteristics, as defined in the <a href=" " target="_blank">Pv5 paper</a>. An eighth group "Unassigned" contains samples from countries with < 25 samples which did not clearly group with another subpopulation. These are colour-coded as follows:
""")
        st.markdown("""
<ul style="list-style-type:none;">
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#4daf4a; border-radius:50%;"></span> LA - Latin America</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#e31a1c; border-radius:50%;"></span> AF-ETH - Africa (Ethiopia)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#984ea3; border-radius:50%;"></span> AS-W - Asia (West)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#9ecae1; border-radius:50%;"></span> AS-SE-W - South-East Asia (West)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#3182bd; border-radius:50%;"></span> AS-SE-E - South-East Asia (East)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#02818a; border-radius:50%;"></span> AS-SE-M - South-East Asia (Maritime)</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#f781bf; border-radius:50%;"></span> OC-NG - Oceania, New Guinea</li>
    <li><span style="display:inline-block; width:10px; height:10px; background-color:#D3D3D3; border-radius:50%;"></span> Unassigned</li>
</ul>
""", unsafe_allow_html=True)

        _st_justify_markdown_html("""
On the x-axis of the geographic distribution subplot of the Haplotype UpSet plot, you will also see the names of lab strains which are also of that haplotype (e.g., P01). 

Due to the country-level aggregation used in the world map plot, countries containing more than one sub-population were allocated their majority sub-population (i.e. Thailand). 
""")
        st.divider()

        st.markdown("## See what's new")
        present_changelog()
        st.divider()

        _st_justify_markdown_html("""
## How to cite

When publishing work that uses data and/or plots from the Pf- or Pv-HaploAtlas, please cite our manuscript on Bioinformatics: 

> Chiyun Lee, Eyyüb S Ünlü, Nina F D White, Jacob Almagro-Garcia, Cristina V Ariani, Richard D Pearson, Pf-HaploAtlas: an interactive web app for spatiotemporal analysis of Plasmodium falciparum genes, Bioinformatics, Volume 40, Issue 11, November 2024, btae673, https://doi.org/10.1093/bioinformatics/btae673
""")
        st.divider()

        _st_justify_markdown_html("""
## Community Quotes""")

        with st.expander("Click to see community quotes"):
            _st_justify_markdown_html("""
"The Pf-HaploAtlas is a legit and nice tool for anyone interested in malaria genomics" from Professor Olivo Miotto, Mahidol Oxford Tropical Medicine Research Unit, Thailand

"We have been using the HaploAtlas app with another PhD student and find it very useful, so thank you to you and the developers of that app! [...] I like being able to download the data to dig in a bit deeper" from Dr Lisa Ranford-Cartwright, University of Glasgow, UK

"When a gene of interest is published in the literature, HaploAtlas is my go-to database for checking its sequence polymorphism. The tool is very simple to use and provides a crucial piece of information" from Associate Professor Antoine Claessens, l'Université de Montpellier, France

"It is an excellent platform that provides valuable information for tracking the evolution and dissemination of P. falciparum drug resistance globally, regionally and locally" from Dr Mary Oboh, Postdoctoral fellow, MRCG at LSHTM, The Gambia

"Pf-HaploAtlas is a powerful yet intuitive tool for exploring global _P. falciparum variation_" from Senior Group Leader Dr Angela Early, Broad Institute, US

""")
        st.divider()

        _st_justify_markdown_html("""
## Acknowledgements

Pv-HaploAtlas currently uses data generated using the MalariaGEN Pv5 data release which was made possible by clinical parasite samples contributed by partner studies, whose investigators are represented in the data release's author list.

""")
        st.divider()

        _st_justify_markdown_html("""
## Contact us
If you'd like to report a bug, request a feature, or give us feedback, please use the following!

- support@malariagen.net
- [our GitHub page](https://github.com/malariagen/pf-haploatlas/issues)
- [this Google Form](https://forms.gle/mDwYr2cPL37dDzPs6)
""")
        st.divider()

        st.markdown("## Created by")
        _show_images_with_urls(
            ["app/files/logo_malariagen.png"],
            ["https://www.malariagen.net/"],
            [70],
            [100]
        )

        st.markdown("## Funded by")
        _show_images_with_urls(
            ["app/files/logo_gates.png"],
            ["https://www.gatesfoundation.org/"],
            [70],
            [100]
        )

        st.divider()

        st.markdown(
            """
<div style="text-align:center">
    Copyright © MalariaGEN at the Liverpool School of Tropical Medicine.
</div>
            """
        , unsafe_allow_html = True)
