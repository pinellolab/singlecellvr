import os
import pathlib
import re

import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
from dash.dependencies import Input, Output, State
import cufflinks as cf

# Initialize app

app = dash.Dash(
    __name__,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1.0"}
    ],
)
server = app.server

# Load data

APP_PATH = str(pathlib.Path(__file__).parent.resolve())

df_lat_lon = pd.read_csv(
    os.path.join(APP_PATH, os.path.join("data", "lat_lon_counties.csv"))
)
df_lat_lon["FIPS "] = df_lat_lon["FIPS "].apply(lambda x: str(x).zfill(5))

df_full_data = pd.read_csv(
    os.path.join(
        APP_PATH, os.path.join("data", "age_adjusted_death_rate_no_quotes.csv")
    )
)
df_full_data["County Code"] = df_full_data["County Code"].apply(
    lambda x: str(x).zfill(5)
)
df_full_data["County"] = (
    df_full_data["Unnamed: 0"] + ", " + df_full_data.County.map(str)
)

YEARS = [2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015]

BINS = [
    "0-2",
    "2.1-4",
    "4.1-6",
    "6.1-8",
    "8.1-10",
    "10.1-12",
    "12.1-14",
    "14.1-16",
    "16.1-18",
    "18.1-20",
    "20.1-22",
    "22.1-24",
    "24.1-26",
    "26.1-28",
    "28.1-30",
    ">30",
]

DEFAULT_COLORSCALE = [
    "#f2fffb",
    "#bbffeb",
    "#98ffe0",
    "#79ffd6",
    "#6df0c8",
    "#69e7c0",
    "#59dab2",
    "#45d0a5",
    "#31c194",
    "#2bb489",
    "#25a27b",
    "#1e906d",
    "#188463",
    "#157658",
    "#11684d",
    "#10523e",
]

DEFAULT_OPACITY = 0.8

mapbox_access_token = "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNqdnBvNDMyaTAxYzkzeW5ubWdpZ2VjbmMifQ.TXcBE-xg9BFdV2ocecc_7g"
mapbox_style = "mapbox://styles/plotlymapbox/cjvprkf3t1kns1cqjxuxmwixz"

# App layout

app.layout = html.Div(
    id="root",
    children=[
        html.Div(
            id="header",
            children=[
                html.Img(id="logo", src=app.get_asset_url("SCVR_logo.png"),style={'height':'15%', 'width':'15%'}),
                html.H4(children="Single Cell VR"),
                html.P(
                    id="description",
                    children="† Single Cell Virtual Reality is a web tool for single-cell data exploration.\
                    VR Headset device required. Google Cardboard is currently the primarily supported headset for this data visualization software.",
                ),
            ],
        ),
        html.Div(
            id="app-container",
            children=[
                html.Div(
                    id="left-column",
                    children=[
                        # html.Div(
                        #     id="slider-container",
                        #     children=[
                        #         html.P(
                        #             id="slider-text",
                        #             children="Drag the slider to change the year:",
                        #         ),
                        #         dcc.Slider(
                        #             id="years-slider",
                        #             min=min(YEARS),
                        #             max=max(YEARS),
                        #             value=min(YEARS),
                        #             marks={
                        #                 str(year): {
                        #                     "label": str(year),
                        #                     "style": {"color": "#7fafdf"},
                        #                 }
                        #                 for year in YEARS
                        #             },
                        #         ),
                        #     ],
                        # ),
                        html.Div(
                            id="dropdown-container",
                            children=[
                                html.P(
                                    id="slider-text",
                                    children="Choose dataset:",
                                ),
                                dcc.Dropdown(
                                    id='chart-dropdown',
                                    options=[
                                        {'label': 'Dataset1', 'value': 'Dataset1'},
                                        {'label': 'Dataset2', 'value': 'Dataset2'},
                                        {'label': 'Dataset3', 'value': 'Dataset3'}
                                    ],
                                    value='Dataset1'
                                ),
                                html.Div(id='dd-output-container'),
                            ],
                        ),
                        html.Div(
                            id="heatmap-container2",
                            children=[
                                html.P("Upload your data:",
                                    id="heatmap-title"),
                                dcc.Upload(
                                    id='upload-data',
                                    children=html.Div([
                                        'Drag and Drop or ',
                                        html.A('Select Files')
                                    ]),
                                    style={
                                        'width': '100%',
                                        'height': '60px',
                                        'lineHeight': '60px',
                                        'borderWidth': '1px',
                                        'borderStyle': 'dashed',
                                        'borderRadius': '5px',
                                        'textAlign': 'center',
                                        'margin': '5px'
                                    },
                                    # Allow multiple files to be uploaded
                                    multiple=True
                                ),
                                html.Div(id='output-data-upload'),
                            ],
                        )
                        # html.Div(
                        #     id="heatmap-container",
                        #     children=[
                        #         html.P(
                        #             "Heatmap of age adjusted mortality rates \
                        #     from poisonings in year {0}".format(
                        #                 min(YEARS)
                        #             ),
                        #             id="heatmap-title",
                        #         ),
                        #         dcc.Graph(
                        #             id="county-choropleth",
                        #             figure=dict(
                        #                 data=[
                        #                     dict(
                        #                         lat=df_lat_lon["Latitude "],
                        #                         lon=df_lat_lon["Longitude"],
                        #                         text=df_lat_lon["Hover"],
                        #                         type="scattermapbox",
                        #                     )
                        #                 ],
                        #                 layout=dict(
                        #                     mapbox=dict(
                        #                         layers=[],
                        #                         accesstoken=mapbox_access_token,
                        #                         style=mapbox_style,
                        #                         center=dict(
                        #                             lat=38.72490, lon=-95.61446
                        #                         ),
                        #                         pitch=0,
                        #                         zoom=3.5,
                        #                     ),
                        #                     autosize=True,
                        #                 ),
                        #             ),
                        #         ),
                        #     ],
                        # ),
                    ],
                ),
                html.Div(
                    id="right-column",
                    children=[
                        html.Div(
                            id="graph-container2",
                            children=[
                                html.P(id="chart-selector", children="Enter VR World:"),
                                html.Div([
                                    # html.Div(dcc.Input(id='input-box', type='text')),
                                    html.A(html.Button("Let's fly!", id='button'),href='http://singlecellvr.com/'),
                                ]),
                            ],
                        ),
                        html.Div(
                            id="video-container",
                            children=[
                                html.P(id="chart-selector2", children="Video tutorial:"),
                                # html.Div([
                                #     # html.Div(dcc.Input(id='input-box', type='text')),
                                #     html.A(html.Button("Let's fly!", id='button'),href='http://singlecellvr.com/'),
                                # ]),
                                html.Div([
                                    # html.Div(dcc.Input(id='input-box', type='text')),
                                    html.Iframe(id='demo_video',src="https://www.youtube.com/embed/i4Zt3JZejbg", width="640",height="380"),
                                ]),
                            ],
                        ),
                    ],
                ),
            ],
        ),
    ],
)


@app.callback(
    dash.dependencies.Output('dd-output-container', 'children'),
    [dash.dependencies.Input('chart-dropdown', 'value')])
def update_output(value):
    return 'You have selected "{}"'.format(value)

# @app.callback(
#     Output("county-choropleth", "figure"),
#     [Input("years-slider", "value")],
#     [State("county-choropleth", "figure")],
# )
def display_map(year, figure):
    cm = dict(zip(BINS, DEFAULT_COLORSCALE))

    data = [
        dict(
            lat=df_lat_lon["Latitude "],
            lon=df_lat_lon["Longitude"],
            text=df_lat_lon["Hover"],
            type="scattermapbox",
            hoverinfo="text",
            marker=dict(size=5, color="white", opacity=0),
        )
    ]

    annotations = [
        dict(
            showarrow=False,
            align="right",
            text="<b>Age-adjusted death rate<br>per county per year</b>",
            font=dict(color="#2cfec1"),
            bgcolor="#1f2630",
            x=0.95,
            y=0.95,
        )
    ]

    for i, bin in enumerate(reversed(BINS)):
        color = cm[bin]
        annotations.append(
            dict(
                arrowcolor=color,
                text=bin,
                x=0.95,
                y=0.85 - (i / 20),
                ax=-60,
                ay=0,
                arrowwidth=5,
                arrowhead=0,
                bgcolor="#1f2630",
                font=dict(color="#2cfec1"),
            )
        )

    if "layout" in figure:
        lat = figure["layout"]["mapbox"]["center"]["lat"]
        lon = figure["layout"]["mapbox"]["center"]["lon"]
        zoom = figure["layout"]["mapbox"]["zoom"]
    else:
        lat = (38.72490,)
        lon = (-95.61446,)
        zoom = 3.5

    layout = dict(
        mapbox=dict(
            layers=[],
            accesstoken=mapbox_access_token,
            style=mapbox_style,
            center=dict(lat=lat, lon=lon),
            zoom=zoom,
        ),
        hovermode="closest",
        margin=dict(r=0, l=0, t=0, b=0),
        annotations=annotations,
        dragmode="lasso",
    )

    base_url = "https://raw.githubusercontent.com/jackparmer/mapbox-counties/master/"
    for bin in BINS:
        geo_layer = dict(
            sourcetype="geojson",
            source=base_url + str(year) + "/" + bin + ".geojson",
            type="fill",
            color=cm[bin],
            opacity=DEFAULT_OPACITY,
            # CHANGE THIS
            fill=dict(outlinecolor="#afafaf"),
        )
        layout["mapbox"]["layers"].append(geo_layer)

    fig = dict(data=data, layout=layout)
    return fig


# @app.callback(Output("heatmap-title", "children"), [Input("years-slider", "value")])
# def update_map_title(year):
#     return "Heatmap of age adjusted mortality rates \
# 				from poisonings in year {0}".format(
#         year
#     )


# @app.callback(
#     Output("selected-data", "figure"),
#     [
#         # Input("county-choropleth", "selectedData"),
#         Input("chart-dropdown", "value"),
#         Input("years-slider", "value"),
#     ],
# )
def display_selected_data(selectedData, chart_dropdown, year):
    if selectedData is None:
        return dict(
            data=[dict(x=0, y=0)],
            layout=dict(
                title="Click-drag on the map to select counties",
                paper_bgcolor="#1f2630",
                plot_bgcolor="#1f2630",
                font=dict(color="#2cfec1"),
                margin=dict(t=75, r=50, b=100, l=75),
            ),
        )
    pts = selectedData["points"]
    fips = [str(pt["text"].split("<br>")[-1]) for pt in pts]
    for i in range(len(fips)):
        if len(fips[i]) == 4:
            fips[i] = "0" + fips[i]
    dff = df_full_data[df_full_data["County Code"].isin(fips)]
    dff = dff.sort_values("Year")

    regex_pat = re.compile(r"Unreliable", flags=re.IGNORECASE)
    dff["Age Adjusted Rate"] = dff["Age Adjusted Rate"].replace(regex_pat, 0)

    if chart_dropdown != "death_rate_all_time":
        title = "Absolute deaths per county, <b>1999-2016</b>"
        AGGREGATE_BY = "Deaths"
        if "show_absolute_deaths_single_year" == chart_dropdown:
            dff = dff[dff.Year == year]
            title = "Absolute deaths per county, <b>{0}</b>".format(year)
        elif "show_death_rate_single_year" == chart_dropdown:
            dff = dff[dff.Year == year]
            title = "Age-adjusted death rate per county, <b>{0}</b>".format(year)
            AGGREGATE_BY = "Age Adjusted Rate"

        dff[AGGREGATE_BY] = pd.to_numeric(dff[AGGREGATE_BY], errors="coerce")
        deaths_or_rate_by_fips = dff.groupby("County")[AGGREGATE_BY].sum()
        deaths_or_rate_by_fips = deaths_or_rate_by_fips.sort_values()
        # Only look at non-zero rows:
        deaths_or_rate_by_fips = deaths_or_rate_by_fips[deaths_or_rate_by_fips > 0]
        fig = deaths_or_rate_by_fips.iplot(
            kind="bar", y=AGGREGATE_BY, title=title, asFigure=True
        )

        fig_layout = fig["layout"]
        fig_data = fig["data"]

        fig_data[0]["text"] = deaths_or_rate_by_fips.values.tolist()
        fig_data[0]["marker"]["color"] = "#2cfec1"
        fig_data[0]["marker"]["opacity"] = 1
        fig_data[0]["marker"]["line"]["width"] = 0
        fig_data[0]["textposition"] = "outside"
        fig_layout["paper_bgcolor"] = "#1f2630"
        fig_layout["plot_bgcolor"] = "#1f2630"
        fig_layout["font"]["color"] = "#2cfec1"
        fig_layout["title"]["font"]["color"] = "#2cfec1"
        fig_layout["xaxis"]["tickfont"]["color"] = "#2cfec1"
        fig_layout["yaxis"]["tickfont"]["color"] = "#2cfec1"
        fig_layout["xaxis"]["gridcolor"] = "#5b5b5b"
        fig_layout["yaxis"]["gridcolor"] = "#5b5b5b"
        fig_layout["margin"]["t"] = 75
        fig_layout["margin"]["r"] = 50
        fig_layout["margin"]["b"] = 100
        fig_layout["margin"]["l"] = 50

        return fig

    fig = dff.iplot(
        kind="area",
        x="Year",
        y="Age Adjusted Rate",
        text="County",
        categories="County",
        colors=[
            "#1b9e77",
            "#d95f02",
            "#7570b3",
            "#e7298a",
            "#66a61e",
            "#e6ab02",
            "#a6761d",
            "#666666",
            "#1b9e77",
        ],
        vline=[year],
        asFigure=True,
    )

    for i, trace in enumerate(fig["data"]):
        trace["mode"] = "lines+markers"
        trace["marker"]["size"] = 4
        trace["marker"]["line"]["width"] = 1
        trace["type"] = "scatter"
        for prop in trace:
            fig["data"][i][prop] = trace[prop]

    # Only show first 500 lines
    fig["data"] = fig["data"][0:500]

    fig_layout = fig["layout"]

    # See plot.ly/python/reference
    fig_layout["yaxis"]["title"] = "Age-adjusted death rate per county per year"
    fig_layout["xaxis"]["title"] = ""
    fig_layout["yaxis"]["fixedrange"] = True
    fig_layout["xaxis"]["fixedrange"] = False
    fig_layout["hovermode"] = "closest"
    fig_layout["title"] = "<b>{0}</b> counties selected".format(len(fips))
    fig_layout["legend"] = dict(orientation="v")
    fig_layout["autosize"] = True
    fig_layout["paper_bgcolor"] = "#1f2630"
    fig_layout["plot_bgcolor"] = "#1f2630"
    fig_layout["font"]["color"] = "#2cfec1"
    fig_layout["xaxis"]["tickfont"]["color"] = "#2cfec1"
    fig_layout["yaxis"]["tickfont"]["color"] = "#2cfec1"
    fig_layout["xaxis"]["gridcolor"] = "#5b5b5b"
    fig_layout["yaxis"]["gridcolor"] = "#5b5b5b"

    if len(fips) > 500:
        fig["layout"][
            "title"
        ] = "Age-adjusted death rate per county per year <br>(only 1st 500 shown)"

    return fig


if __name__ == "__main__":
    app.run_server(debug=True)
