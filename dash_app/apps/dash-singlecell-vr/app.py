import os
import pathlib
import re

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
from dash.dependencies import Input, Output, State
from flask import Flask, send_from_directory,redirect,render_template
from urllib.parse import quote as urlquote
import base64
import uuid
  
# Load data

APP_PATH = str(pathlib.Path(__file__).parent.resolve())

# DATASET_DIRECTORY = os.path.join(APP_PATH, "app_datasets")
UPLOAD_DIRECTORY = os.path.join(APP_PATH, "app_uploaded_files")


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
    external_stylesheets=[dbc.themes.BOOTSTRAP]
)
server = app.server

# @app.before_request
# def before_request():
#     if request.url.startswith('https://'):
#         url = request.url.replace('https://', 'http://', 1)
#         code = 301
#         return redirect(url, code=code)

@server.route("/download/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)

@app.server.route('/view/<uid>')
def serve_static(uid):
	return render_template('index.html', name=uid)

@app.server.route('/help/')
def show_help():
    return render_template('help.html')

# @app.server.route('/view/<uid>')
# def serve_static2(uid):
# 	return render_template('index.html', name=uid)
	
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

app.layout = dbc.Container(
    id="root",
    fluid=True, 
    children=[
        dbc.Row(
            id="header", className="justify-content-between align-items-center flex-nowrap",
            children=[
                html.Img(id='logo', className="col-xs-6 col-sm-6 col-lg-4", src=app.get_asset_url("SCVR_logo.png")),
                html.A(className="col-xs-2 col-sm-2 col-lg-2", href='/help/', children=[dbc.Button("Help", id='button-help', color="dark", disabled=False,n_clicks=0)]),
            ],
        ), 
        dbc.Row(
            children=[
                html.Div( 
                    className="col-container col-lg-4 col-md-12 col-sm-12",
                    children=[
                        html.Div(
                            id="dropdown-container",
                            className="col-content",
                            children=[ 
                                html.H3(
                                    id="slider-text",
                                    children="Choose dataset:",
                                ),
                                dcc.Dropdown(
                                    id='chart-dropdown',
                                    options=[
                                        {'label': 'Mouse blood developmental trajectories', 'value': 'Nestorowa2016-STREAM'},
                                        {'label': 'Mouse myeloid and erythroid differentiation graph', 'value': 'Paul2015-PAGA'},
                                    ],
                                    value=None
                                ),
                                html.Div(id='dd-output-container'),
                                html.Div(id='intermediate-value2', style={'display': 'none'})       
                            ],
                        ),
                    ]),
                html.Div(className="clearfix visible-xs-block visible-sm-block visible-md-block visible-lg-block"),
                dbc.Col(
                    className="col-container",
                    children=[
                        html.Div(
                            id="graph-container2",
                            className="col-content",
                            children=[
                                html.H3(id="chart-selector", children="Enter VR World:"),
                                html.Div([
                                    # html.Div(dcc.Input(id='input-box', type='text')),
                                    # html.Button("Let's fly!", id='button',disabled=False,n_clicks=0),
                                    # html.A(html.Button("Let's fly!", id='button',disabled=False,n_clicks=0),href='/view/123'),
                                    html.Div(id='output-container-button', children=''),
                                    # html.Div(id='output-container-button2',children='')
                                    # html.Div(dcc.Input(id='input-box', type='text')),
                                    # html.A(html.Button("Let's fly!", id='button'),href='http://singlecellvr.com/'),
                                ]),
                            ],
                        ),
                    ]),
            ]),
        dbc.Row( 
            children=[
                dbc.Col(
                    className="col-container col-lg-4 col-sm-12 col-xs-12",
                    children=[
                        html.Div(
                            id="heatmap-container2",
                            className="col-content",
                            children=[
                                html.H3("Or upload your data:",
                                    id="heatmap-title"),
                                dcc.Upload(
                                    id='upload-data',
                                    className="dropper",
                                    children=html.Div(className="dropper-text", children=[
                                        'Drag and Drop or ',
                                        html.A('Select Files')
                                    ]),
                                    # Allow multiple files to be uploaded
                                    multiple=True
                                ),
                                html.Div(id='output-data-upload'),
                                # html.P("Uploaded files:",id="heatmap-title2"),
                                html.Ul(id="file-list"),
                                html.Div(id='intermediate-value', style={'display': 'none'}),
                                html.Div(
                                    children=[
                                        html.P("How to prepare for your submission:",id="heatmap-title2"),
                                        dcc.Markdown('''
                                            * [Generate STREAM trajectories](https://nbviewer.jupyter.org/github/pinellolab/singlecellvr/blob/master/dash_app/apps/dash-singlecell-vr/assets/reformat_files_STREAM.ipynb?flush_cache=true)
                                            * [Generate PAGA graph](https://github.com/pinellolab/singlecellvr/blob/master/dash_app/apps/dash-singlecell-vr/assets/reformat_files_PAGA.ipynb?flush_cache=true)
                                        ''')
                                ])
                        ])]),
                dbc.Col(
                    className="col-container",
                    children=[
                        html.Div(
                            id="video-container",
                            className="col-content col-sm-12 col-xs-12",
                            children=[
                                html.H3(id="chart-selector2", children="Video tutorial:"),
                                html.Div([
                                    html.Iframe(id='demo_video',src="https://www.youtube.com/embed/i4Zt3JZejbg"),
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
    [Output('dd-output-container', 'children'),Output("intermediate-value2", "children")],
    [Input('chart-dropdown', 'value')])
def update_output(value):
    if(value != None):
        if(value=="Nestorowa2016-STREAM"):
            file_id = 'nestorowa2016_stream_report'
        elif(value=="Paul2015-PAGA"):
            file_id = 'paul2015_paga_report'
        return ['You have selected "{}"'.format(value),file_id]
    else:
        return ['No dataset is selected yet',None]

# @app.callback(
#     Output('output-container-button2', 'children'),
#     [Input('intermediate-value2', 'children')])
# def update_output(unique_id):
#     # files = uploaded_files()
#     if unique_id=='':
#         # return 'no files yet'
#         return html.Button("Let's fly!", id='button',disabled=True,n_clicks=0)
#     else:
#         return html.A(html.Button("Let's fly!", id='button',disabled=False,n_clicks=0),href="/view/"+str(unique_id)),

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
    [Input('intermediate-value', 'children'),Input('intermediate-value2', 'children')])
def update_output(unique_id,file_id):
	# files = uploaded_files()

	if(unique_id==None and file_id==None):
		# return 'no files yet'
		return dbc.Button("Please choose or upload dataset!", id='button', className="fly-button", color="link", disabled=True,n_clicks=0)
	if(unique_id !=None):
		return html.A(dbc.Button("Let's fly!", id='button',disabled=False,n_clicks=0, color="link", className="fly-button"),href="/view/"+str(unique_id))
	if(file_id !=None):     
		return html.A(dbc.Button("Let's fly!", id='button',disabled=False,n_clicks=0, color="link", className="fly-button"),href="/view/"+str(file_id))
 
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
    	return [[html.P("")],None]


if __name__ == "__main__":
    # app.run_server(port=os.environ['PORT'])
    app.run_server(port=8050,debug=True,host='0.0.0.0',ssl_context='adhoc')
