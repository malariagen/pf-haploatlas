"""
cd to this directory
run `python3 -m streamlit run app_nina_edits_081123.py` in command line (you may have to pip install streamlit first)
copy the `Network URL` and paste it into your browser (should be something like http://10.160.20.64:8501)
"""

# ============================================================================================================================================================
# Imports 
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

# ============================================================================================================================================================
## Set up title and intros 
st.title('Mutation Discovery App')
st.subheader("Haplotype summaries for *Plasmodium falciparum* genes across time and space")
st.divider()
st.markdown("**Search for a gene in the dropdown list to get started.**")

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

# Read in plot data files and create a dictionary with text between "prep_plot_data" and ".pkl.xz"
file_directory = {
    f: f.replace("prep_plot_data_", "").replace(".pkl.xz", "")
    for f in os.listdir("../array_job_all")
    if f.endswith(".pkl.xz")}

# Create a list of unique gene names and sort it alphabetically
gene_names = sorted(set(file_directory.values()))

# Create a selectbox with the full list of gene names
selected_gene_name = st.selectbox("Select gene name", options=gene_names)
#selected_gene_name = st.selectbox("Select gene name", options=[name.split('_', 2)[-1] for name in gene_names])


# Get the selected full filename based on the selected gene name
filename = [k for k, v in file_directory.items() if v == selected_gene_name][0]

# Retrieve the data from the files
@st.cache_data
def cache_load(gene_id: str):
    with lzma.open(f'../array_job_all/{filename}', 'rb') as file:
        loaded_plot_data = pickle.load(file)
    return loaded_plot_data

df_haplotypes, df_clonal, df_join, background_ns_changes, gene_name = cache_load(filename)

# Add options to change plot settings
with st.expander("Click to edit plot settings"):
    min_samples = st.number_input("Minimum number of samples per haplotype for analysis", min_value = 1, value = 25)
    sample_count_mode = st.radio("Select a mode for Plot 1:", ["Sample counts", "Sample counts on a log scale"], index=0)

# ============================================================================================================================================================
# Create the plots

# Set up layout 
fig = make_subplots(rows = 3, cols = 1, shared_xaxes = True, row_heights = [2, 3, 3], vertical_spacing = 0.025)

# Inputs for plots
total_samples = df_haplotypes['Total'].sum()
df_haplotypes['cum_proportion'] = df_haplotypes['Total'].cumsum() / total_samples 
df_haplotypes_set = df_haplotypes.loc[df_haplotypes['Total'] >= min_samples] 
#df_haplotypes_set['ns_changes'] = df_haplotypes_set['ns_changes'].replace('', '3D7') # Fill in haplotype with no ns_changes as '3D7'

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

# Callback function to update the y-axis based on the radio button selection in previous section
# and change the behaviour of the traces - only show numbers relevant to the y-axis mode chosen by user
if sample_count_mode == "Sample counts":
    y_axis_data = df_haplotypes_set['Total']
    text_labels = y_axis_data
else:
    y_axis_data = np.log(df_haplotypes_set['Total'])
    text_labels = [f'{val:.2f}' for val in y_axis_data]  # Format to 2 decimal places

# Plot 1 - sample counts per haplotype
fig.add_trace(
    go.Bar(
        x=df_haplotypes_set['ns_changes'],
        y=y_axis_data,
        text=text_labels,
        textposition='auto',
        showlegend=False,
        hoverinfo = 'x'
    )
)


# ============================================================================================================================================================
# ============================================================================================================================================================
# Plot 2 - haplotypes across populations

# Define colours for populations
population_colours = collections.OrderedDict()
population_colours['SA']     = "#4daf4a"
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

    # This if statement is an attempt to only display pops with >0 in the trace - not working thus far
    # Check if the maximum proportion value is greater than 0
    if proportion.max() > 0:
        bars.append(go.Bar(
            x=df_haplotypes_set['ns_changes'].values,
            y=proportion,
            marker=dict(color=population_colours[pop]),
            name=pop,
            hovertemplate='%{y:0.1f}%<extra></extra>'  # Use a custom hover template to display percentages with '%' symbol
        ))
        
fig.add_traces(bars, rows=2, cols=1)
fig.update_layout(barmode='stack', legend=dict(x=1, y=0.7))



# ============================================================================================================================================================
# ============================================================================================================================================================
# Plot 3 - UpSet plot

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
            marker=dict(size=18)
        ),
                       rows = 3, cols = 1)
        fig.add_traces(go.Scatter(
            x = [i] * len(other_mutations),
            y = other_mutations,
            showlegend=False,
            marker=dict(size=18)
        ),
                       rows = 3, cols = 1)
    i = i + 1

fig.update_xaxes(tickmode='array', tickvals=[], range = [-0.5, df_haplotypes_set.index.size - 0.5], zeroline = False,
                 showgrid = False, row = 3, col = 1)

fig.update_yaxes(showgrid = True, zeroline = False, gridcolor='rgba(0, 0, 0, 0.15)', 
                 tickvals=df_mutations_set.reset_index()["index"],
                 ticktext=df_mutations_set.reset_index().mutation,
                 row = 3, col = 1)

# ============================================================================================================================================================
# ============================================================================================================================================================

selection_dict = plotly_events(fig, override_height = 900)

if selection_dict == []:
    st.stop()

ns_changes = selection_dict[0]["x"]
if isinstance(ns_changes, int):
    ns_changes = df_haplotypes_set.ns_changes.values[ns_changes]

# ============================================================================================================================================================
# ============================================================================================================================================================

# locations with at least 25 samples (see 20210218_drm_frequency_over_time_pops_adminlevel1.ipynb) (borderline exceptions are marked on this notebook)
# I think we should be able to generate this with code, but for now this is fine
locations = [

    ('SA', 'Colombia', 'Cauca'),

    ('AF-W', 'Gabon', 'Wouleu-Ntem'),
    ('AF-W', 'Cameroon', 'Sud-Ouest'),
    ('AF-W', 'Nigeria', 'Lagos'),
    ('AF-W', 'Benin', 'Littoral'),
    ('AF-W', 'Benin', 'Atlantique'),
    ('AF-W', 'Ghana', 'Volta'),
    ('AF-W', 'Ghana', 'Greater Accra'),
    ('AF-W', 'Ghana', 'Upper East'),
    ('AF-W', 'Ghana', 'Central'),
    ('AF-W', 'Ghana', 'Ashanti'),
    ('AF-W', 'Ghana', 'Brong Ahafo'),
    ('AF-W', 'Burkina Faso', 'Haut-Bassins'),
    ('AF-W', 'Mali', 'Segou'),
    ('AF-W', 'Mali', 'Sikasso'),
    ('AF-W', 'Mali', 'Koulikoro'),
    ('AF-W', 'Mali', 'Bamako'),
    ('AF-W', 'Mali', 'Kayes'),
    
    
    ('AF-W', "Côte d'Ivoire", 'Abidjan'),
    ('AF-W', 'Mauritania', 'Hodh el Gharbi'),
    ('AF-W', 'Mauritania', 'Hodh ech Chargui'), 
    ('AF-W', 'Guinea', 'Nzerekore'),
    ('AF-W', 'Guinea', 'Faranah'),
    ('AF-W', 'Senegal', 'Sedhiou'),
    ('AF-W', 'Senegal', 'Dakar'),
    ('AF-W', 'Gambia', 'Upper River'),
    ('AF-W', 'Gambia', 'North Bank'),
    ('AF-W', 'Gambia', 'Western'),
    
    ('AF-C', 'DRC', 'Kinshasa'),

    ('AF-NE', 'Kenya', 'Kisumu'),
    ('AF-NE', 'Sudan', 'Khartoum'),

    ('AF-E', 'Kenya', 'Kilifi'),
    ('AF-E', 'Mozambique', 'Gaza'), 
    ('AF-E', 'Tanzania', 'Lindi'),
    ('AF-E', 'Tanzania', 'Tanga'),
    ('AF-E', 'Tanzania', 'Morogoro'),
    ('AF-E', 'Tanzania', 'Kagera'),
    ('AF-E', 'Tanzania', 'Kigoma'),
    ('AF-E', 'Malawi', 'Zomba'),   
    ('AF-E', 'Malawi', 'Chikwawa'),

    ('AS-S-E', 'India', 'West Bengal'),
    ('AS-S-E', 'India', 'Odisha'),

     ('AS-S-FE', 'Bangladesh', 'Chittagong'),
     ('AS-S-FE', 'India', 'Tripura'),

     ('AS-SE-W', 'Thailand', 'Tak'),
    ('AS-SE-W', 'Myanmar', 'Tanintharyi'),
    ('AS-SE-W', 'Myanmar', 'Shan'),
    ('AS-SE-W', 'Myanmar', 'Kayin'),
    ('AS-SE-W', 'Myanmar', 'Kachin'),
    ('AS-SE-W', 'Myanmar', 'Bago'),
    ('AS-SE-W', 'Myanmar', 'Mandalay'),
    ('AS-SE-W', 'Myanmar', 'Sagaing'),

    ('AS-SE-E', 'Vietnam', 'Khanh Hoa'),
    ('AS-SE-E', 'Vietnam', 'Ninh Thuan'),
    ('AS-SE-E', 'Vietnam', 'Gia Lai'),
    ('AS-SE-E', 'Vietnam', 'Dak Lak'),
    ('AS-SE-E', 'Vietnam', 'Quang Nam'),
    ('AS-SE-E', 'Vietnam', 'Dak Nong'),
    ('AS-SE-E', 'Vietnam', 'Quang Tri'),
    ('AS-SE-E', 'Vietnam', 'Binh Phuoc'),
    ('AS-SE-E', 'Cambodia', 'Ratanakiri'),
    ('AS-SE-E', 'Cambodia', 'Stueng Traeng'),
    ('AS-SE-E', 'Cambodia', 'Preah Vihear'),
    ('AS-SE-E', 'Cambodia', 'Pursat'),
    ('AS-SE-E', 'Cambodia', 'Battambang'),
    ('AS-SE-E', 'Cambodia', 'Pailin'),

    ('AS-SE-E', 'Laos', 'Sekong'),
    ('AS-SE-E', 'Laos', 'Attapeu'),
    ('AS-SE-E', 'Laos', 'Salavan'),
    ('AS-SE-E', 'Laos', 'Champasak'),
    ('AS-SE-E', 'Laos', 'Savannakhet'),
    ('AS-SE-E', 'Thailand', 'Sisakhet'),

    ('OC-NG', 'Papua New Guinea', 'Milne Bay'),
    ('OC-NG', 'Papua New Guinea', 'Madang'),
    ('OC-NG', 'Papua New Guinea', 'East Sepik'),
    ('OC-NG', 'Indonesia', 'Papua')
]       



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


# df_samples_with_ns_changes = df_samples.join(haplotype_calls['df_all_haplotypes']['ns_changes'])
# Df_join contains the same info as df_samples and haplotype_calls[df_all_haplotypes]
columns_to_remove = ['aa_haplotype', 'nucleotide_haplotype']
df_samples_with_ns_changes = df_join.drop(columns=columns_to_remove, axis=1)
df_samples_with_ns_changes = df_samples_with_ns_changes.loc[df_samples_with_ns_changes['QC pass']]
df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Democratic Republic of the Congo', ['Country']] = 'DRC'

df_samples_with_ns_changes['ns_changes_homozygous'] = ( df_samples_with_ns_changes['ns_changes'] == df_samples_with_ns_changes['ns_changes'].str.upper() )

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

# fig = make_subplots(rows = 3, cols = 1, shared_xaxes = True, row_heights = [2, 3, 2], vertical_spacing = 0)
fig = make_subplots(rows = 1, cols = 3, shared_yaxes = True, column_widths = [1, 1, 5])

xlims = [(1982, 1986),
        (1994,1998),
        (2000,2019)]

for i in [1, 2, 3]:
    fig.update_xaxes(range=xlims[i-1], row=1, col=i)
    
    for pop in reversed(populations):
        
        country = df_frequencies.loc[df_frequencies['Population'] == pop, 'Country']
        _, idx = np.unique(country, return_index=True)
        country = country.reset_index(drop=True)

        for c in country[np.sort(idx)]:

            data = df_frequencies[(df_frequencies['Population'] == pop) & (df_frequencies['Country'] == c)]

            for index, row in data.iterrows():
                if row[ns_changes + ' frequency'] == 0:
                    fig.add_traces(go.Scatter(
                        x = [row.Year],
                        y = [row.Label],
                        showlegend = False,
                        marker=dict(color = "white", size = 18, symbol = "circle-x",
                                    line=dict(color='black',width=2))
                    ), rows = 1, cols = i)
                    
                elif row[ns_changes + ' frequency'] == 1:
                    fig.add_traces(go.Scatter(
                        x = [row.Year],
                        y = [row.Label],
                        showlegend = False,
                        marker=dict(color = "black", size = 18, symbol = "circle-x",
                                    line=dict(color='white',width=2))
                    ), rows = 1, cols = i)
                elif np.isnan(row[ns_changes + ' frequency']):
                    fig.add_traces(go.Scatter(
                        x = [row.Year],
                        y = [row.Label],
                        showlegend = False,
                        marker=dict(color = "white", size = 18, symbol = "circle",
                                    line=dict(color='white',width=2))
                    ), rows = 1, cols = i)
                else:
                    fig.add_traces(go.Scatter(
                        x = [row.Year],
                        y = [row.Label],
                        showlegend = False,
                        marker=dict(color = "white", size = 18, symbol = "circle",
                                    line=dict(color='black',width=2))
                    ), rows = 1, cols = i)
                    fig.add_traces(go.Scatter(
                        x = [row.Year],
                        y = [row.Label],
                        showlegend = False,
                        marker=dict(color = "black", size = 18, symbol = "circle",
                                    opacity=row[ns_changes + ' frequency'])
                    ), rows = 1, cols = i)

fig.update_layout(height = 1100,
                  title = f"Viewing haplotype: {ns_changes}",
                  title_font_size = 18
                  # yaxis=dict(showticklabels=False)
                 )


st.plotly_chart(fig)
