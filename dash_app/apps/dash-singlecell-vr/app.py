import os
import pathlib
import re

import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
from dash.dependencies import Input, Output, State
import cufflinks as cf
from flask import Flask, send_from_directory,redirect,render_template
from urllib.parse import quote as urlquote
import base64
import uuid

# Load data

APP_PATH = str(pathlib.Path(__file__).parent.resolve())

UPLOAD_DIRECTORY = os.path.join(APP_PATH, "app_uploaded_files")
AFRAME_DIRECTORY = os.path.join(APP_PATH,"templates")

# "./dash_app/apps/dash-singlecell-vr/app_uploaded_files"

if not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)

# Normally, Dash creates its own Flask server internally. By creating our own,
# we can create a route for downloading files directly:
# server = Flask(__name__)
# app = dash.Dash(server=server)


# Initialize app

app = dash.Dash(
    __name__,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1.0"}
    ],
)
server = app.server

@server.route("/download/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)

@app.server.route('/view/<uid>')
def serve_static(uid):
	return render_template('index.html', name=uid)
# df_lat_lon = pd.read_csv(
#     os.path.join(APP_PATH, os.path.join("data", "lat_lon_counties.csv"))
# )
# df_lat_lon["FIPS "] = df_lat_lon["FIPS "].apply(lambda x: str(x).zfill(5))

# df_full_data = pd.read_csv(
#     os.path.join(
#         APP_PATH, os.path.join("data", "age_adjusted_death_rate_no_quotes.csv")
#     )
# )
# df_full_data["County Code"] = df_full_data["County Code"].apply(
#     lambda x: str(x).zfill(5)
# )
# df_full_data["County"] = (
#     df_full_data["Unnamed: 0"] + ", " + df_full_data.County.map(str)
# )


# mapbox_access_token = "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNqdnBvNDMyaTAxYzkzeW5ubWdpZ2VjbmMifQ.TXcBE-xg9BFdV2ocecc_7g"
# mapbox_style = "mapbox://styles/plotlymapbox/cjvprkf3t1kns1cqjxuxmwixz"

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
                                html.P("Or upload your data:",
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
                                # html.P("Uploaded files:",id="heatmap-title2"),
                                html.Ul(id="file-list"),
                                # html.Div(id='intermediate-value')
                                html.Div(id='intermediate-value', style={'display': 'none'})
                            ],
                        )
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
                                    # html.Button("Let's fly!", id='button',disabled=False,n_clicks=0),
                                    # html.A(html.Button("Let's fly!", id='button',disabled=False,n_clicks=0),href='/view/123'),
                                    html.Div(id='output-container-button',children='')
                                    # html.Div(dcc.Input(id='input-box', type='text')),
                                    # html.A(html.Button("Let's fly!", id='button'),href='http://singlecellvr.com/'),
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

def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    unique_id = str(uuid.uuid1())
    with open(os.path.join(UPLOAD_DIRECTORY, unique_id+'.zip'), "wb") as fp:
        fp.write(base64.decodebytes(data))
    return unique_id

def uploaded_files():
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(UPLOAD_DIRECTORY):
        path = os.path.join(UPLOAD_DIRECTORY, filename)
        if ((os.path.isfile(path)) and (not filename.startswith('.'))):
            files.append(filename)
    return files

def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "/download/{}".format(urlquote(filename))
    return html.A(filename, href=location)


@app.callback(
    dash.dependencies.Output('dd-output-container', 'children'),
    [dash.dependencies.Input('chart-dropdown', 'value')])
def update_output(value):
    return 'You have selected "{}"'.format(value)


# @app.callback(
#     dash.dependencies.Output('output-container-button', 'children'),
#     [dash.dependencies.Input('button', 'n_clicks')])
# def update_output(n_clicks, value):
#     return 'The input value was "{}" and the button has been clicked {} times'.format(
#         value,
#         n_clicks
#     )


# @app.callback(
#     [Output('button','disabled'),Output('output-container-button', 'children')],
#     [Input('button', 'n_clicks'),Input('intermediate-value', 'children')])
# def update_output(n_clicks,unique_id):
# 	files = uploaded_files()
# 	if len(files) == 0 or n_clicks<2:
# 		return [False,'no files yet']
# 	else:
# 		# return [False,[html.Li(file_download_link(filename+unique_id)) for filename in files]]
# 		return redirect('/aframe/index.html')

# @app.callback(
#     Output('output-container-button', 'children'),
#     [Input('button', 'n_clicks'),Input('intermediate-value', 'children')])
# def update_output(n_clicks,unique_id):
# 	# files = uploaded_files()
# 	if n_clicks<1:
# 		return 'no files yet'
# 	else:
# 		return html.A(unique_id, href="/aframe/index.html?view={}".format(urlquote(unique_id)))
	# else:
	# 	# return [False,[html.Li(file_download_link(filename+unique_id)) for filename in files]]
	# 	return redirect('/aframe/index.html')

@app.callback(
    Output('output-container-button', 'children'),
    [Input('intermediate-value', 'children')])
def update_output(unique_id):
	# files = uploaded_files()
	if unique_id=='':
		# return 'no files yet'
		return html.Button("Let's fly!", id='button',disabled=True,n_clicks=0)
	else:
		return html.A(html.Button("Let's fly!", id='button',disabled=False,n_clicks=0),href="/view/"+str(unique_id)),
		# html.A(unique_id, href="/view/"+str(unique_id))

@app.callback(
    [Output("file-list", "children"),Output("intermediate-value", "children")],
    [Input("upload-data", "filename"), Input("upload-data", "contents")],
)
def update_output(uploaded_filenames, uploaded_file_contents):
    """Save uploaded files and regenerate the file list."""

    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, data in zip(uploaded_filenames, uploaded_file_contents):
            unique_id = save_file(name, data)
            return [[html.P('File ' + name + ' has been uploaded.')],unique_id]
    else:
    	return [[html.P("No files yet!")],'']
    # files = uploaded_files()
    # if len(files) == 0:
    #     return [html.Li("No files yet!")]
    # else:
    #     return [html.Li(file_download_link(filename)) for filename in files]


if __name__ == "__main__":
    app.run_server(host='https://singlecellvr.herokuapp.com', port=os.environ['PORT'], debug=True)
