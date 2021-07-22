import os
import matplotlib as mpl
import pathlib
import re

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
from dash.dependencies import Input, Output, State
from flask import Flask, send_from_directory,redirect,render_template, jsonify, request
import json
from urllib.parse import quote as urlquote
import base64
import uuid
import qrcode
from glob import glob
import requests


APP_PATH = str(pathlib.Path(__file__).parent.resolve())

UPLOAD_DIRECTORY = os.path.join(APP_PATH, "app_uploaded_files")
QR_DIRECTORY = os.path.join(APP_PATH, "assets")
API = 'http://singlecellvrbackend.herokuapp.com'
API = 'http://172.17.0.1:8999'

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
        {"name": "viewport", "content": "width=device-width, initial-scale=1, shrink-to-fit=no, maximum-scale=1.0, user-scalable=no"}
    ],
    external_stylesheets=[dbc.themes.BOOTSTRAP]
)
server = app.server


@server.route("/download/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)


@app.server.route('/view')
def serve_static():
	return render_template('index.html') #, name=uid)


@app.server.route('/help/')
def show_help():
    return render_template('help.html')


app.title = "SingleCellVR"
app.layout = dbc.Container(
    id="root",
    fluid=True,
    children=[
        dbc.Row(
            id="header", className="justify-content-between align-items-center",
            children=[
                html.Img(id='logo', className="col-lg-4 col-md-4 col-sm-12", src=app.get_asset_url("SCVR_logo.png")),
                html.A(className="col-lg-2 col-md-2 col-sm-2", href='/help/', children=[dbc.Button("Help", id='button-help', color="dark", disabled=False, n_clicks=0)]),
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
                                    options=[],
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
                                    html.Div(id='output-container-button', children=''),
                                ]),
                            ],
                        ),
                    ]),
            ]),
        dbc.Row(
            children=[
                html.Div(
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
                                        html.H3("How to prepare for your submission:",id="heatmap-title2"),
                                        dcc.Markdown('''
                                            Check out our package [scvr](https://pypi.org/project/scvr/)
                                            ''')
                                    ]),
                                html.Div(id="qr-output"),
                                # Evil hack
                                html.Div(id="qr-output-secondary")
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
                                    html.Iframe(id='demo_video',src="https://www.youtube.com/embed/CcEpnrw34zM"),
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
    Output('chart-dropdown', 'options'),
    [Input('dropdown-container', 'n_clicks')]
)
def update_options(n_clicks):
    options = json.loads(requests.get(f'{API}/databases').text)
    return options


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
    Output(component_id='qr-output', component_property='children'),
    [Input(component_id='intermediate-value', component_property='children')]
)
def render_qrcode(unique_id):
    if unique_id:
        save_qr_image(unique_id)
        return html.Div(
                    children=[
                        html.H3("Uploaded dataset:"),
                        html.Img(src='/assets/' + str(unique_id) + '.bmp', style={'width': '40%', 'margin-bottom': '20px'}),
                    ]
                )

def save_qr_image(unique_id):
    img = qrcode.make("https://singlecellvr.herokuapp.com/view/" + str(unique_id))
    i = img.get_image()
    if unique_id:
        i.save(os.path.join(QR_DIRECTORY) + '/' + str(unique_id) + '.bmp')

@app.callback(
    [Output('dd-output-container', 'children'),
     Output('intermediate-value2', 'children'), 
     Output('qr-output-secondary', 'children')],
    [Input('chart-dropdown', 'value')])
def update_output(value):
    if(value != None):
        file_id = ''
        if(value=="Nestorowa16"):
            file_id = 'Nestorowa16_Stream'
        elif(value=="Paul2015"):
            file_id = 'Paul2015_Paga'
        elif(value=="Macosko2015"):
            file_id = 'Macosko2015_Scanpy'
        elif(value=="10xPBMC"):
            file_id = '10xPBMC_Seurat'
        elif(value=="Pancrease"):
            file_id = 'Pancrease_Velocity'
        elif(value=='nestorowa'):
            file_id = 'nestorowa_Velocity'
        save_qr_image(file_id)
        print("Value " + str(value))
        return ['You have selected "{}"'.format(value),
                file_id, 
                html.Div(
                    children=[
                        html.H3("Preprocessed dataset:"),
                        html.Img(src='/assets/' + file_id + '.bmp', style={'width': '40%'}),
                    ]
                )]
    else:
        return ['No dataset is selected yet',None,None]

@app.callback(
    Output('output-container-button', 'children'),
    [Input('intermediate-value', 'children'),Input('intermediate-value2', 'children')])
def update_output(unique_id,file_id):
    if not unique_id and not file_id:
        return dbc.Button("Please choose or upload dataset!", id='button', className="fly-button", color="link", disabled=True,n_clicks=0)
    if unique_id:
        return html.A(dbc.Button("Let's fly!", id='button',disabled=False,n_clicks=0, color="link", className="fly-button"),
                                    href="/view?dataset="+str(unique_id) + "&fulldataset=false")
    if file_id:     
        return html.A(dbc.Button("Let's fly!", id='button',disabled=False,n_clicks=0, color="link", className="fly-button"),
                                    href="/view?dataset="+str(file_id) + "&fulldataset=true")
  
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

