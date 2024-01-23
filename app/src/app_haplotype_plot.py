import streamlit as st

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from plotly.subplots import make_subplots
from streamlit_plotly_events import plotly_events

from src.utils import cache_load_population_colours

def generate_haplotype_plot(df_haplotypes, gene_id_selected, background_ns_changes, min_samples, sample_count_mode):
    
    population_colours = cache_load_population_colours()
    
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
    
    return ns_changes, df_haplotypes_set