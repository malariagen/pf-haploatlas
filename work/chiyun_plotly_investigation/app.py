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
# ============================================================================================================================================================

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

file_directory = {f: f.replace("prep_plot_data_", "").replace(".pkl.xz", "") for f in os.listdir("../produce_gene_summaries") if f.endswith(".pkl.xz")}
filename = st.selectbox("Select file", options = file_directory)
gene_id = file_directory[filename]

@st.cache_data
def cache_load(gene_id: str):
    with lzma.open(f'../produce_gene_summaries/prep_plot_data_{gene_id}.pkl.xz', 'rb') as file:
        loaded_plot_data = pickle.load(file)
        
    return loaded_plot_data

df_haplotypes, df_clonal, df_join, background_ns_changes, haplotype_calls, gene_name, df_samples = cache_load(gene_id)

with st.expander("Click to edit settings"):
    cumulative_mode = st.checkbox("Cumulative statistics", value = False)
    log_y_scale = st.checkbox("Y-axis log scale", value = False)
    min_samples = st.number_input("Minimum number of samples per haplotype for analysis", min_value = 1, value = 25)
    

# ============================================================================================================================================================
# ============================================================================================================================================================


fig = make_subplots(rows = 3, cols = 1, shared_xaxes = True, row_heights = [2, 3, 2], vertical_spacing = 0)

total_samples = df_haplotypes['Total'].sum()
df_haplotypes['cum_proportion'] = df_haplotypes['Total'].cumsum() / total_samples
df_haplotypes_set = df_haplotypes.loc[df_haplotypes['Total'] >= min_samples]

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

fig.add_trace(
    go.Bar(
        x = df_haplotypes_set['ns_changes'],
        y = df_haplotypes_set['Total'],
        text = df_haplotypes_set['Total'],
        textposition = 'auto',
        showlegend=False
    ),
    row = 1, col = 1
)

# ============================================================================================================================================================
# ============================================================================================================================================================

bars = []
for column_num in np.arange(2, df_haplotypes_set.shape[1]-4):
    pop = df_haplotypes_set.columns[column_num]
    
    bars += [go.Bar(x = df_haplotypes_set.ns_changes.values,
                    y = df_haplotypes_set[df_haplotypes_set.columns[column_num]] / df_haplotypes_set.Total,
                    marker = dict(color = population_colours[pop]),
                    name = pop,
                   )]
    
fig.add_traces(bars, rows = 2, cols = 1)

fig.update_layout(barmode = 'stack', legend=dict(x=1, y=0.7))

# ============================================================================================================================================================
# ============================================================================================================================================================


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
            showlegend=False
            ,
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

fig.update_yaxes(showgrid = True, zeroline = False,
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

#     ('SA', 'Colombia', 'Cauca'),

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
    
#     ('AS-S-E', 'India', 'West Bengal'),
#     ('AS-S-E', 'India', 'Odisha'),

#     ('AS-S-FE', 'Bangladesh', 'Chittagong'),
#     ('AS-S-FE', 'India', 'Tripura'),

#     ('AS-SE-W', 'Thailand', 'Tak'),
#     ('AS-SE-W', 'Myanmar', 'Tanintharyi'),
#     ('AS-SE-W', 'Myanmar', 'Shan'),
#     ('AS-SE-W', 'Myanmar', 'Kayin'),
#     ('AS-SE-W', 'Myanmar', 'Kachin'),
#     ('AS-SE-W', 'Myanmar', 'Bago'),
#     ('AS-SE-W', 'Myanmar', 'Mandalay'),
#     ('AS-SE-W', 'Myanmar', 'Sagaing'),

#     ('AS-SE-E', 'Vietnam', 'Khanh Hoa'),
#     ('AS-SE-E', 'Vietnam', 'Ninh Thuan'),
#     ('AS-SE-E', 'Vietnam', 'Gia Lai'),
#     ('AS-SE-E', 'Vietnam', 'Dak Lak'),
#     ('AS-SE-E', 'Vietnam', 'Quang Nam'),
#     ('AS-SE-E', 'Vietnam', 'Dak Nong'),
#     ('AS-SE-E', 'Vietnam', 'Quang Tri'),
#     ('AS-SE-E', 'Vietnam', 'Binh Phuoc'),
#     ('AS-SE-E', 'Cambodia', 'Ratanakiri'),
#     ('AS-SE-E', 'Cambodia', 'Stueng Traeng'),
#     ('AS-SE-E', 'Cambodia', 'Preah Vihear'),
#     ('AS-SE-E', 'Cambodia', 'Pursat'),
#     ('AS-SE-E', 'Cambodia', 'Battambang'),
#     ('AS-SE-E', 'Cambodia', 'Pailin'),
    
#     ('AS-SE-E', 'Laos', 'Sekong'),
#     ('AS-SE-E', 'Laos', 'Attapeu'),
#     ('AS-SE-E', 'Laos', 'Salavan'),
#     ('AS-SE-E', 'Laos', 'Champasak'),
#     ('AS-SE-E', 'Laos', 'Savannakhet'),
#     ('AS-SE-E', 'Thailand', 'Sisakhet'),

#     ('OC-NG', 'Papua New Guinea', 'Milne Bay'),
#     ('OC-NG', 'Papua New Guinea', 'Madang'),
#     ('OC-NG', 'Papua New Guinea', 'East Sepik'),
#     ('OC-NG', 'Indonesia', 'Papua')
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


df_samples_with_ns_changes = df_samples.join(haplotype_calls['df_all_haplotypes']['ns_changes'])
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
                  # yaxis=dict(showticklabels=False)
                 )

st.dataframe(df_frequencies)
    
st.plotly_chart(fig)