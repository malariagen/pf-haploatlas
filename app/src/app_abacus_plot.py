import streamlit as st
import collections
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from src.utils import cache_load_population_colours, generate_download_buttons, _cache_load_utility_mappers, _st_justify_markdown_html

def _plotly_arrow(x0, x1, y):
    """One time function used to generate the arrow in the legend of the abacus plot"""
    a = go.layout.Annotation(
        x = x1, ax = x0, y = y, ay = y,
        xref="x", yref="y", text="", showarrow=True,
        axref="x", ayref='y', arrowhead=3, arrowwidth=1.5)

    return a

def _locations_agg(x, ns_changes):
    """Aggregation function used to reformat dataframe in preparation for abacus plot"""
    names = collections.OrderedDict()
    names['n'] = np.count_nonzero(x['ns_changes_homozygous'])
    if names['n'] == 0:
        names[f'{ns_changes} frequency'] = np.nan
    else:
        names[f'{ns_changes} frequency'] = np.count_nonzero(
            ( x['ns_changes'] == ns_changes)
        ) / names['n']

    return pd.Series(names)

def _partial_frequency_marker_colour(freq: float) -> str:
    """
    Convenience function which takes haplotype frequency and returns an rgba string for the grey colour used to
    colour the marker in the abacus plot. The higher the frequency, the darker the grey
    """
    colour_intensity = max(0, min(255, 255 - int(255 * freq)))
    marker_colour = f"rgba({colour_intensity}, {colour_intensity}, {colour_intensity}, 1)"

    return marker_colour


def generate_abacus_plot(ns_changes, df_join, min_samples, df_haplotypes_set, gene_id_selected):
    """Main function called in main.py to generate and present the abacus plot"""
    
    population_colours = cache_load_population_colours()
    utility_mappers = _cache_load_utility_mappers()

    gene_name_selected = utility_mappers["gene_ids_to_gene_names"][gene_id_selected]

    # Filter QC fail and missing samples
    df_samples_with_ns_changes = df_join.loc[
        (df_join['QC pass'] == True) &
        (df_join["Exclusion reason"].isin(['Analysis_set', 'Unassigned_population'])) &
        (df_join["HaploAtlas exclusion reason"].isin(['Analysis_set', 'Unassigned_population']))
    ].copy()
    
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Democratic Republic of the Congo', ['Country']] = 'DRC'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes.ns_changes == "", "ns_changes"] = "P01"

    df_samples_with_ns_changes["Admin level 1"] = df_samples_with_ns_changes["Admin level 1"].fillna("(Unspecified)")

    df_samples_with_ns_changes.loc[:,'ns_changes_homozygous'] = ( df_samples_with_ns_changes['ns_changes'] == df_samples_with_ns_changes['ns_changes'].str.upper() )

    aggregated_locations = df_samples_with_ns_changes.groupby(['Population', 'Country', 'Admin level 1']).apply(lambda x: len(x) >= min_samples)
    locations = aggregated_locations.index[aggregated_locations.values].values

    df_frequencies = (
        df_samples_with_ns_changes
        .groupby(['Population', 'Year', 'Country', 'Admin level 1'])
        .apply(lambda x: _locations_agg(x, ns_changes))
        .reset_index()
        .set_index(['Population', 'Country', 'Admin level 1'])
        .loc[locations]
        .reset_index()
    )
    df_frequencies = df_frequencies.loc[df_frequencies['n'] >= min_samples]
    df_frequencies['Label'] = df_frequencies['Country'] + ', ' + df_frequencies['Admin level 1']

    populations = [col for col in df_haplotypes_set.columns if col not in ["number_of_mutations", "ns_changes", "Total", "ns_changes_list", "sample_names", "cum_proportion"]]

    # ============================================================================================================================================================
    # ============================================================================================================================================================

    st.divider()

    st.subheader(f'2. Abacus plot: {ns_changes}')
    _st_justify_markdown_html("""
The Abacus plot shows how the haplotype frequency of your selected haplotype changes across locations and time (in years). The colour intensity of each “bead” on the Abacus plot corresponds to its observed frequency in each year and in each location. Haplotypes at fixation are marked with “100” to highlight 100% frequency, whilst “beads” with 0% haplotype frequency are crossed out. Hover your mouse over the data to see details. 

Click and drag to zoom to focus on certain locations. Double-click to reset. 
""")
    
    fig = make_subplots(rows = 2, cols = 2,
                        vertical_spacing = 0,
                        horizontal_spacing = 0,
                        specs = [
                            [{"colspan": 2}, None],
                            [{}, {}],
                        ],
                        row_heights = [1, 20],
                        column_widths = [1, 3],
                        shared_yaxes = "rows"
                       )

    labels_list = []
    population_colours_list = []

    def _abacus_scatter(**kwargs):
        """Convenience function for creating scatter points on the abacus plot"""
        return go.Scatter(
            customdata = [[row.n,
                           np.rint(row.n * row[ns_changes + ' frequency']), # int casting occasionally rounded floats to 0.0, so using `np.rint`
                           np.round(row[ns_changes + ' frequency'] * 100, 1)]],
            showlegend = False,
            **kwargs
        )

    scatter_config = {
        "zero_frequency": {
            "marker": dict(color="white",
                           size=11,
                           symbol="circle-x",
                           line=dict(color='gray',
                                     width=1))
        },
        "full_frequency": {
            "marker": dict(color="black",
                           size=16,
                           symbol="circle"),
            "text": "<b>100</b>",
            "mode": "markers+text",
            "textfont": dict(size=8,
                             color="white")
        },
    }

    hovertemplate = '<b>%{y} in %{x}</b><br>Samples with selected haplotype: %{customdata[1]} (%{customdata[2]}%)<br>Number of samples: %{customdata[0]}<extra></extra>',
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
                                        ), rows = 2, cols = 2)

                elif row[ns_changes + ' frequency'] == 1:
                    fig.add_traces(
                        _abacus_scatter(x = [row.Year], y = [row.Label],
                                        hovertemplate = hovertemplate,
                                        **scatter_config["full_frequency"]
                                        ), rows = 2, cols = 2)

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
                                        ), rows = 2, cols = 2)
    
    min_year, max_year = df_frequencies.Year.min(), df_frequencies.Year.max()
    n_locations = df_frequencies.Label.nunique()

    fig.update_xaxes(title_text=  "Year", range = (min_year - 1, max_year + 1), row = 2, col = 2)
    fig.update_yaxes(title_text = "Location", range = (-0.5, n_locations - 0.5), row = 2, col = 1)

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
        _abacus_scatter(x = [0.6], y = [legend_y], hoverinfo = "none", **scatter_config["full_frequency"]),
        go.Scatter(x = [0.4, 0.6], y = [0.6, 0.6], hoverinfo = "none", showlegend = False, mode = "text", text = ["0%", "100%"]),
        go.Scatter(x = [0.5], y = [0.9], hoverinfo = "none", showlegend = False, mode = "text", text = ["Haplotype Frequency"])] +

        [_abacus_scatter(x = [pos], y = [legend_y], hoverinfo = "none",
                         marker=dict(color = _partial_frequency_marker_colour(freq),
                                     size = 16, symbol = "circle",
                                     line=dict(color='black', width=1.5))
                        ) for pos, freq in zip(partial_frequency_positions, partial_frequency_frequencies)],
        
        rows = 1, cols = 1)
    
    def map_value(x):
        """
        Linear mapping of number of locations to plot height. I roughly used these coords to fit:
        (17, 1000), (100, 1800), (160, 2500)
        """
        y = 1000 + 10.5 * (x - 17)
        return max(1000, min(2500, y))

    plot_height = map_value(len(df_frequencies[["Country", "Admin level 1"]].drop_duplicates()))

    fig.add_annotation(_plotly_arrow(0.42, 0.57, legend_y+0.3))

    fig.update_layout(
        title={
            'text': f"Pf-HaploAtlas Abacus plot: {gene_name_selected} ({ns_changes})",
            'y':0.99,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {
                'size': 14,
            }},
        height = plot_height,
        xaxis = dict(tickvals = [], range = (0, 1), fixedrange=True, zeroline=False),
        xaxis2 = dict(range = (0, 1), fixedrange=True, tickvals = []),
        xaxis3 = dict(fixedrange=True, tickangle=-60, tickvals = np.arange(min_year, max_year).astype(int)),

        yaxis = dict(tickvals = [], range = (0, 1), fixedrange=True, zeroline=False),
        yaxis2 = dict(fixedrange=True, tickvals = []),
        yaxis3 = dict(fixedrange=True, showticklabels = False, tickmode='linear'),
        margin=dict(t=50, b=55, l=0, r=0)
    )

    st.plotly_chart(fig, config = {"displayModeBar": False})

    generate_download_buttons(fig, gene_id_selected, 1300, 800, plot_number = 2)
