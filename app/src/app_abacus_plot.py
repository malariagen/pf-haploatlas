import streamlit as st
import collections

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from plotly.subplots import make_subplots

from src.utils import cache_load_population_colours

def generate_abacus_plot(ns_changes, df_join, min_samples, df_haplotypes_set):
    population_colours = cache_load_population_colours()
    
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