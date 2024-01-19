# ============================================================================================================================================================
import json
import os
import lzma
import pickle
import collections

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from plotly.subplots import make_subplots
from streamlit_plotly_events import plotly_events
from PIL import Image
from webbrowser import open_new_tab

# ============================================================================================================================================================
## Set up title and intros 
st.set_page_config(
    page_title = "Pf7 mutation discovery app",
    page_icon = Image.open("logo.png"), #PLACEHOLDER
    initial_sidebar_state = "expanded",
    menu_items={
        'Get Help': 'https://plasmodb.org/plasmo/app/record/gene/PF3D7_0523000#PubMed', #PLACEHOLDER
        'Report a bug': "https://plasmodb.org/plasmo/app/record/gene/PF3D7_0523000#PubMed", #PLACEHOLDER
        'About': "PLACEHOLDER"
    }
)

st.title('Pf7.0 Haplotype Explorer (PfHEx7.0)')
st.subheader("Haplotype summaries for *Plasmodium falciparum* genes across time and space")
st.divider()
placeholder = st.empty()
placeholder.markdown("### Search for a gene below to get started.")

# ============================================================================================================================================================
## Sidebar for additional info
with st.sidebar:
    st.title("**Further Information**")
    st.divider()
    st.header("**Samples**")
    st.markdown("The Mutation Discovery App uses 16,203 QC pass samples from the [Pf7 dataset.](https://wellcomeopenresearch.org/articles/8-22/v1)")
    st.header("**Genes**")
    st.markdown("The Mutation Discovery App uses **XXX** genes, according to the **XXX** version of the 3D7 Pf genome from *PlasmoDB link?* All genes have a unique identifier, e.g. **PF3D7_XXXXXXX**, and in some cases a gene name, e.g. **MDR1**.")
    st.header("**Plots**")
    st.markdown("The app generates four plots per gene:")
    st.markdown("**1. Sample counts per haplotype** - cumulative counts of samples per haplotype. Toggle the y-axis between raw values or a log scale.")
    st.markdown("**2. Population proportions per haplotype** - relative contributions of each of the major Pf7 sub-populations to each haplotype.")
    st.markdown("**3. UpSet plot** - view the mutation makeup of each haplotype.")
    st.markdown("**4. Abacus plot** - click on a haplotype to view its frequency across first-level administrative divisions and years with at least 25 samples.")
    st.markdown("The shade of the point represents the haplotype frequency from white (0%) to black (100%). Where frequency is exactly 0% or 100% the point is marked with a cross to represent fixation.")

# ============================================================================================================================================================
# Read in the plot data and provide options to customise the settings

@st.cache_data
def cache_load_gene_mapper():
    with open("gene_mapping.json", "r") as f:
        return json.load(f)

def _is_core_genome(filename: str):
    if "VAR" in filename:
        return False
    elif "SURF" in filename:
        return False
    elif "RIF" in filename:
        return False
    else:
        return True

def _lookup_gene_name(gene_id: str) -> str:
    if gene_id == "--":
        return "--"
    else:
        return gene_ids_to_gene_names[gene_id]

gene_mapper = cache_load_gene_mapper()
default_value = "--"
base_path = "../backend/23-11-23_run_array_job"


files_to_gene_ids = {f: f.split(".")[0] for f in os.listdir(base_path) if f.endswith("pkl.xz")}
gene_ids_to_files = dict(zip(files_to_gene_ids.values(), files_to_gene_ids.keys()))

gene_ids_to_gene_names = {gene: (f'{gene} - {gene_mapper[gene]}' if not gene_mapper[gene] is "." else gene) for gene in sorted(files_to_gene_ids.values())}
gene_ids = [gene_id for gene_id, gene_name in gene_ids_to_gene_names.items() if _is_core_genome(gene_name)]


# Handle URL - functionalise this
current_query_params = st.experimental_get_query_params()
gene_id_extracted = (current_query_params["gene_id"][0] if "gene_id" in current_query_params.keys() else "")

if gene_id_extracted and "gene_id" not in st.session_state:
    st.session_state["gene_id"] = gene_id_extracted

gene_id_selected = st.selectbox("", ["--"] + gene_ids, format_func = _lookup_gene_name, key = "gene_id")

if gene_id_selected == "--":
    placeholder.markdown("### Search for a gene below to get started.")
    st.experimental_set_query_params()
    st.stop()

if gene_id_selected is not gene_id_extracted:
    current_query_params["gene_id"] = gene_id_selected
    st.experimental_set_query_params(**current_query_params)

filename = gene_ids_to_files[gene_id_selected]
placeholder.empty()

# ============================================================================================================================================================
# Retrieve the data from the files
@st.cache_data
def cache_load(gene_id: str):
    with lzma.open(f'{base_path}/{filename}', 'rb') as file:
        loaded_plot_data = pickle.load(file)
    return loaded_plot_data

df_haplotypes, df_clonal, df_join, background_ns_changes, gene_name = cache_load(filename)

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

# ============================================================================================================================================================

# Inputs for plots
total_samples = df_haplotypes['Total'].sum()
df_haplotypes['cum_proportion'] = df_haplotypes['Total'].cumsum() / total_samples 
df_haplotypes_set = df_haplotypes.loc[df_haplotypes['Total'] >= min_samples] 
df_haplotypes_set['ns_changes'] = df_haplotypes_set['ns_changes'].replace('', '3D7 REF')

if len(df_haplotypes_set) == 0:
    st.warning("No haplotype data found.")
    st.stop()

mutations = pd.Series(np.unique(np.concatenate(df_haplotypes_set['ns_changes_list'].values)))
mutations = mutations.loc[mutations != '']
mutations_np_ref = mutations.apply(lambda x: x[1:])
aa = mutations_np_ref.apply(lambda x: int(x[:-1]))
df_mutations_set = pd.DataFrame(
    {
        'mutation': mutations.values,
        'aa': aa.values,
    }
).sort_values('aa').reset_index(drop=True).reset_index().set_index('mutation')

# Some arbitrary plot-scaling calculations
upset_plot_height = 1.5 + len(df_mutations_set) / 5
total_plot_height = (5 + upset_plot_height) * 100

# Create the plots
fig = make_subplots(rows = 3, cols = 1, shared_xaxes = True, row_heights = [2, 3, upset_plot_height], vertical_spacing = 0)

# Plot 1 - sample counts per haplotype
fig.add_trace(
    go.Bar(
        x=df_haplotypes_set['ns_changes'],
        y=df_haplotypes_set['Total'],
        text=df_haplotypes_set['Total'],
        textposition='auto',
        showlegend=False,
        hoverinfo = 'x',
        hovertemplate = "<b>%{x}:</b> %{y}<extra></extra>"
    )
)

if sample_count_mode == "Sample counts on a log scale":
    fig.update_yaxes(type = "log", row = 1, col = 1)


# ============================================================================================================================================================
# ============================================================================================================================================================
# Plot 2 - haplotypes across populations

# Define colours for populations
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

bars = []
for column_num in np.arange(2, df_haplotypes_set.shape[1] - 4):
    pop = df_haplotypes_set.columns[column_num]
    proportion = ((df_haplotypes_set[df_haplotypes_set.columns[column_num]] / df_haplotypes_set['Total']).round(3))*100
    
    bars.append(go.Bar(
        x=df_haplotypes_set['ns_changes'].values,
        y=proportion,
        marker=dict(color=population_colours[pop]),
        customdata = [pop] * len(df_haplotypes_set),
        name=pop,
        hovertemplate='<b>%{customdata}:</b> %{y:0.1f}%<extra></extra>'
    ))

fig.add_traces(bars, rows=2, cols=1)
fig.update_layout(barmode='stack', legend=dict(x=1, y = 1 - 2 / (5 + upset_plot_height)))
# 2 is the row height ratio of the first plot, so this calculation places the legend right below the first plot

# ============================================================================================================================================================
# ============================================================================================================================================================
# Plot 3 - UpSet plot

marker_size = 5 + np.sqrt(len(df_haplotypes_set))

i = 0
if ( not '' in background_ns_changes ) and background_ns_changes in df_haplotypes_set['ns_changes'].values:
    background_mutation_indices = df_mutations_set.loc[
        df_haplotypes_set.loc[
            df_haplotypes_set['ns_changes'] == background_ns_changes,
            'ns_changes_list'
        ].values[0],
        'index'
    ].values
else:
    background_mutation_indices = []

for ix, row in df_haplotypes_set.iterrows():
    if row['ns_changes_list'][0] != '':
        indexes = df_mutations_set.loc[row['ns_changes_list'], 'index'].values
        fig.add_traces(go.Scatter(
            x = [i] * len(indexes),
            y = indexes,
            showlegend=False,
            hoverinfo = 'none',
            mode = "lines"),
                       rows = 3, cols = 1)
        
        
        background_mutations = np.intersect1d(indexes, background_mutation_indices)
        other_mutations = np.setdiff1d(indexes, background_mutation_indices)
        
        fig.add_traces(go.Scatter(
            x = [i] * len(background_mutations),
            y = background_mutations,
            showlegend=False,
            hovertemplate='%{y}<extra></extra>',
            marker=dict(size=marker_size)
        ),
                       rows = 3, cols = 1)
        fig.add_traces(go.Scatter(
            x = [i] * len(other_mutations),
            y = other_mutations,
            showlegend=False,
            hovertemplate='%{y}<extra></extra>',
            marker=dict(size=marker_size)
        ),
                       rows = 3, cols = 1)
    i = i + 1

fig.update_xaxes(tickmode='array', tickvals=[], range = [-0.5, df_haplotypes_set.index.size - 0.5], zeroline = False,
                 showgrid = False, row = 3, col = 1)

fig.update_yaxes(showgrid = True, zeroline = False, gridcolor='rgba(0, 0, 0, 0.15)', 
                 tickvals=df_mutations_set.reset_index()["index"],
                 ticktext=df_mutations_set.reset_index().mutation,
                 row = 3, col = 1)
fig.update_layout(hovermode = 'closest')

# ============================================================================================================================================================
# ============================================================================================================================================================

selection_dict = plotly_events(fig, override_height = total_plot_height)

if selection_dict == []:
    st.stop()

ns_changes = selection_dict[0]["x"]
if isinstance(ns_changes, int):
    ns_changes = df_haplotypes_set.ns_changes.values[ns_changes]

# ============================================================================================================================================================
# ============================================================================================================================================================

def locations_agg(x, ns_changes):
    names = collections.OrderedDict()
    names['n'] = np.count_nonzero(x['ns_changes_homozygous'])
    if names['n'] == 0:
        names[f'{ns_changes} frequency'] = np.nan
    else:
        names[f'{ns_changes} frequency'] = np.count_nonzero(
            ( x['ns_changes'] == ns_changes)
        ) / names['n']
    
    return pd.Series(names)

columns_to_remove = ['aa_haplotype', 'nucleotide_haplotype']
df_samples_with_ns_changes = df_join.drop(columns=columns_to_remove, axis=1)
df_samples_with_ns_changes = df_samples_with_ns_changes.loc[df_samples_with_ns_changes['QC pass']]
df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Democratic Republic of the Congo', ['Country']] = 'DRC'
df_samples_with_ns_changes.loc[df_samples_with_ns_changes.ns_changes == "", "ns_changes"] = "3D7 REF"

df_samples_with_ns_changes['ns_changes_homozygous'] = ( df_samples_with_ns_changes['ns_changes'] == df_samples_with_ns_changes['ns_changes'].str.upper() )

aggregated_locations = df_samples_with_ns_changes.groupby(['Population', 'Country', 'Admin level 1']).apply(lambda x: len(x) >= min_samples)
locations = aggregated_locations.index[aggregated_locations.values].values

df_frequencies = (
    df_samples_with_ns_changes
    .groupby(['Population', 'Year', 'Country', 'Admin level 1'])
    .apply(lambda x: locations_agg(x, ns_changes))
    .reset_index()
    .set_index(['Population', 'Country', 'Admin level 1'])
    .loc[locations]
    .reset_index()
)
df_frequencies = df_frequencies.loc[df_frequencies['n'] >= min_samples]
df_frequencies['Label'] = df_frequencies['Country'] + ', ' + df_frequencies['Admin level 1']

populations = [col for col in df_haplotypes_set.columns if col not in ["number_of_mutations", "ns_changes", "Total", "ns_changes_list", "sample_names", "cum_proportion"]]
label = {p: df_frequencies.loc[df_frequencies['Population'] == p, 'Admin level 1'].nunique() for p in populations}

# ============================================================================================================================================================
# ============================================================================================================================================================


st.divider()

# def _haplotype_formatter(h: str) -> str:
#     if h.count("/") > 5:
#         h = "/".join(h.split("/")[:5]) + "<br>" + "/".join(h.split("/")[5:])
#     else:
#         return h

st.subheader(f"Viewing haplotype: {ns_changes}")
st.write("Click and drag to zoom. Double-click to reset.")
    
fig = make_subplots(rows = 2, cols = 4,
                    vertical_spacing = 0,
                    horizontal_spacing = 0.05,
                    specs = [
                        [{"colspan": 4}, None, None, None],
                        [{}, {}, {}, {}],
                    ],
                    row_heights = [1, 20],
                    column_widths = [3, 1, 1, 5],
                    shared_yaxes = "rows"
                   )


xlims = [
    (1982, 1986),
    (1994,1998),
    (2000,2019)
]

labels_list = []
population_colours_list = []


def _abacus_scatter(**kwargs):
    return go.Scatter(
        customdata = [[row.n, int(row.n * row[ns_changes + ' frequency']), np.round(row[ns_changes + ' frequency'] * 100, 1)]],
        showlegend = False,
        **kwargs
    )

scatter_config = {
    "zero_frequency": {
        "marker": dict(color="white",
                       size=13,
                       symbol="circle-x",
                       line=dict(color='gray',
                                 width=1))
    },
    "full_frequency": {
        "marker": dict(color="black",
                       size=20,
                       symbol="circle"),
        "text": "<b>!</b>",
        "mode": "markers+text",
        "textfont": dict(size=30,
                         color="white")
    },
}

def _partial_frequency_marker_colour(freq: float):
    colour_intensity = max(0, min(255, 255 - int(255 * freq)))
    marker_colour = f"rgba({colour_intensity}, {colour_intensity}, {colour_intensity}, 1)"
    
    return marker_colour
    

hovertemplate = '<b>%{y} in %{x}</b><br>Samples with selected haplotype: %{customdata[1]} (%{customdata[2]}%)<br>All samples: %{customdata[0]}<extra></extra>',

for i in [2, 3, 4]:
    fig.update_xaxes(range=xlims[i-2], row=2, col=i)
    
    for pop in reversed(populations):
        
        country = df_frequencies.loc[df_frequencies['Population'] == pop, 'Country']
        _, idx = np.unique(country, return_index=True)
        country = country.reset_index(drop=True)

        for c in country[np.sort(idx)]:

            data = df_frequencies[(df_frequencies['Population'] == pop) & (df_frequencies['Country'] == c)]
                
            for index, row in data.iterrows():
                
                if row.Label not in labels_list:
                    labels_list.append(row.Label)
                    population_colours_list.append(population_colours[row.Population])
                
                if row[ns_changes + ' frequency'] == 0:
                    fig.add_traces(
                        _abacus_scatter(x = [row.Year], y = [row.Label],
                                        hovertemplate = hovertemplate,
                                        **scatter_config["zero_frequency"]
                                       ), rows = 2, cols = i)
                    
                elif row[ns_changes + ' frequency'] == 1:
                    fig.add_traces(
                        _abacus_scatter(x = [row.Year], y = [row.Label],
                                        hovertemplate = hovertemplate,
                                        **scatter_config["full_frequency"]
                                       ), rows = 2, cols = i)

                else:

                    
                    fig.add_traces(
                        _abacus_scatter(x = [row.Year], y = [row.Label],
                                        hovertemplate = hovertemplate,
                                        marker=dict(color = _partial_frequency_marker_colour(row[ns_changes + ' frequency']),
                                                    size = 16, symbol = "circle",
                                                    line=dict(
                                                        color='black',
                                                        width=1.5)
                                                   )
                                       ), rows = 2, cols = i)
                            
fig.add_traces(
    go.Scatter(
        x = np.ones(len(labels_list)),
        y = labels_list,
        mode="text",
        text = labels_list,
        textposition = "middle left",
        textfont=dict(color=population_colours_list),
        showlegend = False,
        hoverinfo = "none"
    ),
    rows = 2, cols = 1
)

legend_y = 0.3
partial_frequency_frequencies = np.linspace(0.05, 0.95, 8)
partial_frequency_positions = np.linspace(0.45, 0.55, 8)

fig.add_traces([
    _abacus_scatter(x = [0.4], y = [legend_y], hoverinfo = "none", **scatter_config["zero_frequency"]),
    _abacus_scatter(x = [0.6], y = [legend_y], hoverinfo = "none", **scatter_config["full_frequency"])] +
    
    [_abacus_scatter(x = [pos], y = [legend_y], hoverinfo = "none", marker=dict(color = _partial_frequency_marker_colour(freq), size = 16, symbol = "circle",
                                                                                line=dict(color='black', width=1.5))) for pos, freq in zip(partial_frequency_positions, partial_frequency_frequencies)] +
    
    [go.Scatter(x = [0.4, 0.6], y = [0.6, 0.6], hoverinfo = "none", showlegend = False, mode = "text", text = ["0%", "100%"])] +
    
    [go.Scatter(x = [0.5], y = [0.9], hoverinfo = "none", showlegend = False, mode = "text", text = ["Haplotype Frequency"])],
    
    rows = 1, cols = 1)


def plotly_arrow(x0, x1, y):
    a = go.layout.Annotation(
        x = x1, ax = x0, y = y, ay = y,
        xref="x", yref="y", text="", showarrow=True,
        axref="x", ayref='y', arrowhead=3, arrowwidth=1.5)
    
    return a

fig.add_annotation(plotly_arrow(0.42, 0.57, legend_y+0.3))

fig.update_layout(height = 1300,
                  xaxis = dict(tickvals = [], range = (0, 1), fixedrange=True, zeroline=False),
                  xaxis2 = dict(range = (0, 1), fixedrange=True, tickvals = []),
                  xaxis3 = dict(fixedrange=True, tickangle=-60),
                  xaxis4 = dict(fixedrange=True, tickangle=-60),
                  xaxis5 = dict(fixedrange=True, tickangle=-60, tickvals = np.arange(2000, 2020).astype(int)),
                  
                  yaxis = dict(tickvals = [], range = (0, 1), fixedrange=True, zeroline=False),
                  yaxis2 = dict(tickvals = []),
                  yaxis3 = dict(showticklabels = False, tickmode='linear'),
                  yaxis4 = dict(showticklabels = False, tickmode='linear'),
                  yaxis5 = dict(showticklabels = False, tickmode='linear'),
                  
                 )

st.plotly_chart(fig, config = {"displayModeBar": False})