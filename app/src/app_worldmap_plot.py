import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import collections

from src.utils import cache_load_population_colours, generate_download_buttons, _cache_load_utility_mappers

def _locations_agg(x, ns_changes):
    """Aggregation function used to reformat dataframe in preparation for worldmap plot"""
    names = collections.OrderedDict()
    names['n'] = np.count_nonzero(x['ns_changes_homozygous'])
    if names['n'] == 0:
        names[f'frequency'] = 0
        names['haplo_count'] = 0
    else:
        names['haplo_count'] = np.count_nonzero(( x['ns_changes'] == ns_changes))
        names[f'frequency'] = names['haplo_count'] / names['n']

    return pd.Series(names)

def _partial_frequency_marker_colour(freq: float) -> str:
    """
    Convenience function which takes haplotype frequency and returns an rgba string for the grey colour used to
    colour the marker in the worldmap plot. The higher the frequency, the darker the grey
    """
    colour_intensity = max(0, min(255, 255 - int(255 * freq/100)))
    marker_colour = f"rgba({colour_intensity}, {colour_intensity}, {colour_intensity}, 1)"

    return marker_colour

def _abacus_scatter(**kwargs):
        """Convenience function for creating scatter points on the haplotype frequency legend"""
        return go.Scatter(
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
            "text": "<b>100</b>",
            "mode": "markers+text",
            "textfont": dict(size=8,
                             color="white")
        },
    }

def _plotly_arrow(x0, x1, y):
    """One time function used to generate the arrow in the legend of the worldmap plot"""
    a = go.layout.Annotation(
        x = x1, ax = x0, y = y, ay = y,
        xref="x", yref="y", text="", showarrow=True,
        axref="x", ayref='y', arrowhead=3, arrowwidth=1.5)

    return a    

def generate_worldmap_plot(ns_changes, df_join, min_samples, gene_id_selected):
    """Main function called in main.py to generate and present the worldmap plot"""

    st.divider()

    st.subheader(f'3. Worldmap plot: {ns_changes}')
    st.markdown(
        f"""
The Worldmap plot displays the average haplotype frequency over an interval of time (in years) at a country-level on a global map. As above, the colour intensity of each “bead” corresponds to the frequency, and beads are coloured by geographic distribution (see sidebar for details). Hover your mouse over the data to see details. 

Adjust the slider below to choose your time interval of interest for calculating the proportion of samples containing the {ns_changes} haplotype: 
        """
    )
    year = st.slider('', 1982, 2024, (2010, 2018))
    population_colours = cache_load_population_colours()
    utility_mappers = _cache_load_utility_mappers()

    gene_name_selected = utility_mappers["gene_ids_to_gene_names"][gene_id_selected]

    ### Data pre-processing
    # ideally, in the future, this part needs te be handled by data_formatter
    if len(df_join) == 0:
        st.warning("No haplotype data found.")
        st.stop()
    
    # Filter QC fail and missing samples  
    df_join= df_join[df_join['Exclusion reason'] == 'Analysis_set']

    df_samples_with_ns_changes = df_join
    # worldmap map requires iso_alpha values
    plotly_worldmap_df = px.data.gapminder().query("year==2007")
    iso_country_dict = dict(zip(plotly_worldmap_df['country'], plotly_worldmap_df['iso_alpha']))
    df_samples_with_ns_changes['iso_alpha'] = df_samples_with_ns_changes['Country'].map(iso_country_dict)
    df_samples_with_ns_changes['iso_alpha'] = df_samples_with_ns_changes['Country'].map(iso_country_dict)
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Papua New Guinea', 'iso_alpha'] = 'PNG'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Laos', 'iso_alpha'] = 'LAO'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Democratic Republic of the Congo', 'iso_alpha'] = 'COD'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == "Côte d'Ivoire", 'iso_alpha'] = 'CIV'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'South Sudan', 'iso_alpha'] = 'SSD'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == "Lao People's Democratic Republic", 'iso_alpha'] = 'LAO'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'United Republic of Tanzania', 'iso_alpha'] = 'TZA'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'The Gambia', 'iso_alpha'] = 'GMB'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Guyana', 'iso_alpha'] = 'GUY'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'DRC', 'iso_alpha'] = 'COD'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Solomon Islands', 'iso_alpha'] = 'SLB'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Vanuatu', 'iso_alpha'] = 'VUT'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Congo', 'iso_alpha'] = 'COG'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'French Guiana', 'iso_alpha'] = 'GUF'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Yemen', 'iso_alpha'] = 'YEM'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Suriname', 'iso_alpha'] = 'SUR'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Cape Verde', 'iso_alpha'] = 'CPV'
    df_samples_with_ns_changes['iso_alpha'] = df_samples_with_ns_changes['iso_alpha'].astype(object)

    # Fix for population of vietnam
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == "Vietnam", ['Population']] = 'AS-SE-E'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == "India", ['Population']] = 'AF-E'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == "Kenya", ['Population']] = 'AS-S-E'

    # deal with encoding of 'wildtype' in literature dataset
    if 'wildtype' in df_samples_with_ns_changes['ns_changes'].values:
        df_samples_with_ns_changes.loc[df_samples_with_ns_changes['ns_changes'] == 'wildtype', 'ns_changes'] = ''

    df_samples_with_ns_changes = df_samples_with_ns_changes.loc[df_samples_with_ns_changes['QC pass']]
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'Democratic Republic of the Congo', ['Country']] = 'DRC'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == 'United Republic of Tanzania', ['Country']] = 'Tanzania'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes['Country'] == "Lao People's Democratic Republic", ['Country']] = 'Laos'
    df_samples_with_ns_changes.loc[df_samples_with_ns_changes.ns_changes == "", "ns_changes"] = "3D7 REF"

    df_samples_with_ns_changes['ns_changes_homozygous'] = ( df_samples_with_ns_changes['ns_changes'] == df_samples_with_ns_changes['ns_changes'].str.upper() )

    # Use range slider to filter relative years in the dataframe
    df_samples_with_ns_changes= df_samples_with_ns_changes[(df_samples_with_ns_changes['Year']>=year[0]) & (df_samples_with_ns_changes['Year']<=year[1])]
    df_samples_with_ns_changes['Year-interval'] = str(year)

    ### AGGREGATION     
    if len(df_samples_with_ns_changes.groupby(['iso_alpha', 'Country', 'Year-interval', 'Population'])) == 0:
        st.warning("No haplotype data found.")
        st.stop()
        
    df_frequencies = (
        df_samples_with_ns_changes
        .groupby(['iso_alpha', 'Country', 'Year-interval', 'Population'])
        .apply(lambda x: _locations_agg(x, ns_changes))
        .reset_index()
        .set_index(['Country'])
        .reset_index())
    
    # only>min_samples 
    df_frequencies = df_frequencies.loc[(df_frequencies['n'] >= min_samples)]

    df_frequencies['frequency'] = np.round(df_frequencies['frequency']*100,2)
    df_frequencies[['n', 'haplo_count']] = df_frequencies[['n', 'haplo_count']].astype('int')

    ### WORLDMAP PLOT (WORLD MAP)

    # Create a subplot comprising haplotype frequency (legend) and world map 
    fig = make_subplots(rows=2, cols=1, specs=[[{"type": "xy"}], [{"type": "scattergeo"}]],
                        vertical_spacing = 0,
                        row_heights = [1, 10])

    # Add haplotype frequency legend as subplot 1
    legend_y = 0.3  # Adjust this value to position the legend vertically
    partial_frequency_frequencies = np.linspace(5, 95, 8)
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

    fig.add_annotation(_plotly_arrow(0.42, 0.57, legend_y+0.3))
    fig.update_layout(xaxis = dict(tickvals = [], dtick=1, range = (0, 1), 
                                #    fixedrange=True, 
                                   zeroline=False),
                      yaxis = dict(tickvals = [], range = (0, 1),
                                #    fixedrange=True,
                                   zeroline=False))

    # Add worldmap plot (scattergeo subplot)
    for _, row in df_frequencies.iterrows():
        trace = go.Scattergeo(
            locations=[row['iso_alpha']],
            hoverinfo='text',
            hovertemplate=f"<b>{row['Country']}: {row['Year-interval'].strip('()').replace(',', ' - ')}</b><br>Population: {row['Population']}<br>Samples with selected haplotype: {row['haplo_count']} ({row['frequency']}%) <br>Number of samples: {row['n']}</b><extra></extra>",
            marker=dict(
                size=13,
                line=dict(
                    color=population_colours[row['Population']],
                    width=1.35
                )
            ), showlegend=False
        )

        if row['frequency'] == 0:
            trace.marker.symbol = 'circle-x'
            trace.marker.color = 'white'
        elif row['frequency'] == 100:
            trace.marker.symbol = 'circle'
            trace.mode = "markers+text"
            trace.marker.color = 'black'
        else:
            trace.marker.symbol = 'circle'
            trace.marker.color = _partial_frequency_marker_colour(row['frequency'])

        fig.add_trace(trace, row=2, col=1)

    # Update layout
    fig.update_layout(
        title={
            'text': f"Pf-HaploAtlas Worldmap plot: {gene_name_selected} ({ns_changes})",
            'y':0.99,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {
                'size': 14,
            }},
        height=600, width=800,
        margin=dict(t=40, b=5, l=5, r=5)
    )
    fig.update_geos(projection_type="natural earth")

    st.plotly_chart(fig, config = {"displayModeBar": False})

    generate_download_buttons(fig, gene_id_selected, 600, 800, plot_number = 3)