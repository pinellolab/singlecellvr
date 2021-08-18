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
from flask_cors import CORS, cross_origin
import json
from urllib.parse import quote as urlquote
import base64
import uuid
import qrcode
from glob import glob
import requests

import scanpy as sc
from scvr import converters
from scipy.sparse import isspmatrix
from pandas.api.types import is_string_dtype,is_numeric_dtype
import networkx as nx
import numpy as np
import stream as st
import resource
import gc

APP_PATH = str(pathlib.Path(__file__).parent.resolve())
DATASET_DIRECTORY = os.path.join(APP_PATH, "app_datasets")
UPLOAD_DIRECTORY = os.path.join(APP_PATH, "app_uploaded_files")
QR_DIRECTORY = os.path.join(APP_PATH, "assets")
API = 'https://singlecellvr.pinellolab.partners.org'
#API = 'http://0.0.0.0:8080'

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
cors = CORS(server)
server.config['CORS_HEADERS'] = 'Content-Type'

@server.route("/download/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(DATASET_DIRECTORY, path, as_attachment=True)


@app.server.route('/view')
def serve_static():
	return render_template('index.html') #, name=uid)


@app.server.route('/help/')
def show_help():
    return render_template('help.html')

def get_tool_type(file):
    for tool in ['stream', 'paga', 'scanpy', 'seurat']:
        if tool in file.lower():
            return tool

def datasets_payload():
    adata_files = glob(os.path.join(DATASET_DIRECTORY, "*"))
    adata_list = list()
    for file in adata_files:
        tool = get_tool_type(file)
        label = os.path.basename(file).split(".")[0]
        label_list = label.split("_")
        adata_list.append(
            {
                "value": label_list[0],
                "label": label,
                "path": file,
                "type": tool,
            }
        )
    return adata_list

@app.server.route("/databases", methods=["GET", "POST"])
def get_databases():
    return jsonify(dataset_payload)

@app.server.route("/data_type", methods=["GET", "POST"])
def get_dataset_type():
    """
    http://127.0.0.1:8000/data_type?db_name=1_scanpy_10xpbmc
    """
    db_name = request.args.get("db_name")
    return jsonify({"type": db_name.split("_")[1]})


def get_dataset_type_adata(db_name):
    return db_name.split("_")[1]


@app.server.route("/coordinates", methods=["GET", "POST"])
def get_coordinates():
    """
    http://127.0.0.1:8000/coordinates?db_name=1_scanpy_10xpbmc&embed=umap
    http://127.0.0.1:8000/coordinates?db_name=3_velocity_pancrease&embed=umap
    http://127.0.0.1:8000/coordinates?db_name=4_seurat_10xpbmc&embed=umap
    http://127.0.0.1:8000/coordinates?db_name=5_stream_nestorowa16&embed=umap
    """
    db_name = request.args.get("db_name")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0]

    try:
        del adata
    except:
        pass

    if get_dataset_type_adata(db_name).lower() in ["scanpy", "velocity", "seurat", "paga"]:
        adata = sc.read(filename)
        embed = request.args.get("embed")
    else:
        adata = st.read(filename, file_format="pkl", workdir="./")

    list_cells = []
    for i in range(adata.shape[0]):
        dict_coord_cells = dict()
        dict_coord_cells["cell_id"] = adata.obs_names[i]
        if get_dataset_type_adata(db_name).lower() in ["scanpy", "paga", "velocity"]:
            dict_coord_cells["x"] = str(adata.obsm[f"X_{embed}"][i, 0])
            dict_coord_cells["y"] = str(adata.obsm[f"X_{embed}"][i, 1])
            dict_coord_cells["z"] = str(adata.obsm[f"X_{embed}"][i, 2])
        elif get_dataset_type_adata(db_name).lower() == "seurat":
            dict_coord_cells["x"] = str(adata.obsm[f"{embed}_cell_embeddings"][i, 0])
            dict_coord_cells["y"] = str(adata.obsm[f"{embed}_cell_embeddings"][i, 1])
            dict_coord_cells["z"] = str(adata.obsm[f"{embed}_cell_embeddings"][i, 2])
        elif get_dataset_type_adata(db_name).lower() == "stream":
            file_path = os.path.join(adata.uns["workdir"], "test")
            if not os.path.exists(file_path):
                os.makedirs(file_path)
            flat_tree = adata.uns["flat_tree"]
            epg = adata.uns["epg"]
            epg_node_pos = nx.get_node_attributes(epg, "pos")
            ft_node_label = nx.get_node_attributes(flat_tree, "label")
            ft_node_pos = nx.get_node_attributes(flat_tree, "pos")
            list_curves = []
            for edge_i in flat_tree.edges():
                branch_i_pos = np.array(
                    [epg_node_pos[i] for i in flat_tree.edges[edge_i]["nodes"]]
                )
                df_coord_curve_i = pd.DataFrame(branch_i_pos)
                dict_coord_curves = dict()
                dict_coord_curves["branch_id"] = (
                    ft_node_label[edge_i[0]] + "_" + ft_node_label[edge_i[1]]
                )
                dict_coord_curves["xyz"] = [
                    {
                        "x": df_coord_curve_i.iloc[j, 0],
                        "y": df_coord_curve_i.iloc[j, 1],
                        "z": df_coord_curve_i.iloc[j, 2],
                    }
                    for j in range(df_coord_curve_i.shape[0])
                ]
                list_curves.append(dict_coord_curves)

            ## output topology of stream graph
            dict_nodes = dict()
            list_edges = []
            for node_i in flat_tree.nodes():
                dict_nodes_i = dict()
                dict_nodes_i["node_name"] = ft_node_label[node_i]
                dict_nodes_i["xyz"] = {
                    "x": ft_node_pos[node_i][0],
                    "y": ft_node_pos[node_i][1],
                    "z": ft_node_pos[node_i][2],
                }
                dict_nodes[ft_node_label[node_i]] = dict_nodes_i
            for edge_i in flat_tree.edges():
                dict_edges = dict()
                dict_edges["nodes"] = [
                    ft_node_label[edge_i[0]],
                    ft_node_label[edge_i[1]],
                ]
                dict_edges["weight"] = 1
                list_edges.append(dict_edges)

            list_cells = []
            for i in range(adata.shape[0]):
                dict_coord_cells = dict()
                dict_coord_cells['cell_id'] = adata.obs_names[i]
                dict_coord_cells['x'] = adata.obsm['X_dr'][i,0]
                dict_coord_cells['y'] = adata.obsm['X_dr'][i,1]
                dict_coord_cells['z'] = adata.obsm['X_dr'][i,2]
                list_cells.append(dict_coord_cells)
            return jsonify(
                {"nodes": dict_nodes, "edges": list_edges, "graph": list_curves, "cells": list_cells}
            )
        else:
            raise TypeError("not supported format")
        list_cells.append(dict_coord_cells)
    del adata
    gc.collect()
    return jsonify(list_cells)


@app.server.route("/features", methods=["GET", "POST"])
def get_features():
    """
    scanpy examples:
      http://127.0.0.1:8000/features?db_name=1_scanpy_10xpbmc&feature=louvain
      http://127.0.0.1:8000/features?db_name=1_scanpy_10xpbmc&feature=expression&gene=SUMO3

    seurat examples:
      http://127.0.0.1:8000/features?db_name=4_seurat_10xpbmc&feature=expression&gene=SUMO3
      http://127.0.0.1:8000/features?db_name=4_seurat_10xpbmc&feature=expression&gene=SUMO3

    velocity examples:
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=clusters
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=expression&gene=Rbbp7
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=velocity&embed=umap&time=None
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=velocity&embed=umap&time=1
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=velocity&embed=umap&time=10

    velocity grid examples:
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=velocity_grid&embed=umap&time=10
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=velocity_grid&embed=umap&time=100
    """
    database = request.args.get("db_name")
    feature = request.args.get("feature")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{database}.*"))[0]

    db_type = get_dataset_type_adata(filename)
    if feature.lower() == "velocity":
        embed = request.args.get("embed")

    try:
        del adata
    except:
        pass

    if get_dataset_type_adata(database).lower() in ["scanpy", "velocity", "seurat", "paga"]:
        adata = sc.read(filename)
    else:
        adata = st.read(filename, file_format="pkl", workdir="./")

    list_metadata = []
    if feature in get_available_annotations_adata(adata):  # cluster columns
        if f"{feature}_colors" in adata.uns.keys():
            dict_colors = {
                feature: dict(
                    zip(adata.obs[feature].cat.categories, adata.uns[f"{feature}_colors"])
                )
            }
        else:
            dict_colors = {
                feature: dict(
                    zip(adata.obs[feature], converters.get_colors(adata, feature))
                )
            }
        for i in range(adata.shape[0]):
            dict_metadata = dict()
            dict_metadata["cell_id"] = adata.obs_names[i]
            dict_metadata["label"] = adata.obs[feature].tolist()[i]
            dict_metadata["clusters"] = adata.obs[feature].tolist()[i]
            dict_metadata["clusters_color"] = dict_colors[feature][
                dict_metadata["clusters"]
            ]
            list_metadata.append(dict_metadata)
    elif feature in ["expression", "rna"]:  # pseudotime or latent_time columns
        gene = request.args.get("gene")
        if gene not in adata.var_names:
            return jsonify({})
        else:
            if "time" in feature:
                values = adata.obs[feature]
            else:
                if db_type == "seurat":
                    values = (
                        adata[:, gene].layers["norm_data"].toarray()[:, 0]
                        if isspmatrix(adata.layers["norm_data"])
                        else adata[:, gene].layers["norm_data"][:, 0]
                    )
                else:
                    values = (
                        adata[:, gene].X.toarray()[:, 0]
                        if isspmatrix(adata.X)
                        else adata[:, gene].X[:, 0]
                    )

            cm = mpl.cm.get_cmap("viridis", 512)
            norm = mpl.colors.Normalize(vmin=0, vmax=max(values), clip=True)
            list_metadata = []
            for i, x in enumerate(adata.obs_names):
                dict_genes = dict()
                dict_genes["cell_id"] = x
                dict_genes["color"] = mpl.colors.to_hex(cm(norm(values[i])))
                list_metadata.append(dict_genes)
    elif feature == "velocity":
        list_metadata = []
        time = request.args.get("time")
        for i in range(adata.shape[0]):
            dict_coord_cells = dict()
            if isinstance(adata.obs_names[i], bytes):
                dict_coord_cells["cell_id"] = adata.obs_names[i].decode("utf-8")
            else:
                dict_coord_cells["cell_id"] = adata.obs_names[i]

            dict_coord_cells["x"] = str(adata.obsm[f"X_{embed}"][i, 0])
            dict_coord_cells["y"] = str(adata.obsm[f"X_{embed}"][i, 1])
            dict_coord_cells["z"] = str(adata.obsm[f"X_{embed}"][i, 2])

            if time == "None":
                dict_coord_cells["x1"] = str(adata.obsm[f"velocity_{embed}"][i, 0])
                dict_coord_cells["y1"] = str(adata.obsm[f"velocity_{embed}"][i, 1])
                dict_coord_cells["z1"] = str(adata.obsm[f"velocity_{embed}"][i, 2])
            elif time in list(map(str, [0.01, 0.1, 1, 5, 10, 20, 30, 50, 100])):
                dict_coord_cells["x1"] = str(
                    adata.obsm[f"absolute_velocity_{embed}_{time}s"][i, 0]
                )
                dict_coord_cells["y1"] = str(
                    adata.obsm[f"absolute_velocity_{embed}_{time}s"][i, 1]
                )
                dict_coord_cells["z1"] = str(
                    adata.obsm[f"absolute_velocity_{embed}_{time}s"][i, 2]
                )
            else:
                return jsonify({})
            list_metadata.append(dict_coord_cells)
    elif feature == "velocity_grid":
        list_metadata = []
        time = request.args.get("time")
        p_mass = adata.uns['p_mass']
        for i in np.where(p_mass >= 1)[0]:
            dict_coord_cells = dict()

            if time == "None":
                dict_coord_cells["x"] = str(adata.uns[f"X_grid"][i, 0])
                dict_coord_cells["y"] = str(adata.uns[f"X_grid"][i, 1])
                dict_coord_cells["z"] = str(adata.uns[f"X_grid"][i, 2])
                dict_coord_cells["x1"] = str(adata.uns[f"V_grid"][i, 0])
                dict_coord_cells["y1"] = str(adata.uns[f"V_grid"][i, 1])
                dict_coord_cells["z1"] = str(adata.uns[f"V_grid"][i, 2])
            elif time in list(map(str, [0.01, 0.1, 1, 5, 10, 20, 50, 80, 100])):
                dict_coord_cells["x"] = str(adata.uns[f"X_grid_{time}"][i, 0])
                dict_coord_cells["y"] = str(adata.uns[f"X_grid_{time}"][i, 1])
                dict_coord_cells["z"] = str(adata.uns[f"X_grid_{time}"][i, 2])
                dict_coord_cells["x1"] = str(
                    adata.uns[f"V_grid_{time}"][i, 0]
                )
                dict_coord_cells["y1"] = str(
                    adata.uns[f"V_grid_{time}"][i, 1]
                )
                dict_coord_cells["z1"] = str(
                    adata.uns[f"V_grid_{time}"][i, 2]
                )
            else:
                return jsonify({})
            list_metadata.append(dict_coord_cells)
    elif feature == "paga":
        G = nx.from_numpy_matrix(adata.uns["paga"]["connectivities"].toarray())
        adata.uns["paga"]["pos"] = get_paga3d_pos(adata)
        ## output coordinates of paga graph
        list_lines = []
        for edge_i in G.edges():
            dict_coord_lines = dict()
            dict_coord_lines["branch_id"] = [[str(edge_i[0]), str(edge_i[1])]]
            dict_coord_lines["xyz"] = [
                {"x": pos[0], "y": pos[1], "z": pos[2]}
                for pos in adata.uns["paga"]["pos"][[edge_i[0], edge_i[1]], :]
            ]
            list_lines.append(dict_coord_lines)

        ## output topology of paga graph
        dict_nodes = dict()
        list_edges = []
        dict_nodename = {
            i: adata.obs[adata.uns["paga"]["groups"]].cat.categories[i]
            for i in G.nodes()
        }
        for node_i in G.nodes():
            dict_nodes_i = dict()
            dict_nodes_i["node_name"] = dict_nodename[node_i]
            dict_nodes_i["xyz"] = {
                "x": adata.uns["paga"]["pos"][:, 0][node_i],
                "y": adata.uns["paga"]["pos"][:, 1][node_i],
                "z": adata.uns["paga"]["pos"][:, 2][node_i],
            }
            dict_nodes[node_i] = dict_nodes_i
        for edge_i in G.edges():
            dict_edges = dict()
            dict_edges["nodes"] = [str(edge_i[0]), str(edge_i[1])]
            dict_edges["weight"] = adata.uns["paga"]["connectivities"][
                edge_i[0], edge_i[1]
            ]
            list_edges.append(dict_edges)
        list_metadata = {"nodes": dict_nodes, "edges": list_edges}
    elif feature == "curves":
        flat_tree = adata.uns['flat_tree']
        epg = adata.uns['epg']
        epg_node_pos = nx.get_node_attributes(epg,'pos')
        ft_node_label = nx.get_node_attributes(flat_tree,'label')
        ft_node_pos = nx.get_node_attributes(flat_tree,'pos')
        list_curves = []
        for edge_i in flat_tree.edges():
            branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])
            df_coord_curve_i = pd.DataFrame(branch_i_pos)
            dict_coord_curves = dict()
            dict_coord_curves['branch_id'] = ft_node_label[edge_i[0]] + '_' + ft_node_label[edge_i[1]]
            dict_coord_curves['xyz'] = [{'x':df_coord_curve_i.iloc[j,0],
                                         'y':df_coord_curve_i.iloc[j,1],
                                         'z':df_coord_curve_i.iloc[j,2]} for j in range(df_coord_curve_i.shape[0])]
            list_curves.append(dict_coord_curves)
        list_metadata = list_curves
    del adata
    gc.collect()
    return jsonify({feature: list_metadata})

@app.server.route("/legend", methods=["GET", "POST"])
def get_legend_colors():
    database = request.args.get("db_name")
    ann = request.args.get("feature")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{database}.*"))[0]

    db_type = get_dataset_type_adata(filename)
    if ann.lower() == "velocity":
        embed = request.args.get("embed")

    try:
        del adata
    except:
        pass

    if get_dataset_type_adata(database).lower() in ["scanpy", "velocity", "seurat", "paga"]:
        adata = sc.read(filename)
    else:
        adata = st.read(filename, file_format="pkl", workdir="./")

    categories = []
    colors = []
    if(is_numeric_dtype(adata.obs[ann])):
        cm = mpl.cm.get_cmap('viridis',512)
        norm = mpl.colors.Normalize(vmin=0, vmax=max(adata.obs[ann]),clip=True)
        for x in sorted(adata.obs[ann]):
            colors.append(mpl.colors.to_hex(cm(norm(x))))
        categories = adata.obs[ann].tolist()
    else:
        categories = adata.obs[ann].astype('category').cat.categories.tolist()
        length = len(categories)
        # check if default matplotlib palette has enough colors
        # mpl.style.use('default')
        if len(mpl.rcParams['axes.prop_cycle'].by_key()['color']) >= length:
            cc = mpl.rcParams['axes.prop_cycle']()
            palette = [next(cc)['color'] for _ in range(length)]
        else:
            if length <= 20:
                palette = default_20
            elif length <= 28:
                palette = default_28
            elif length <= len(default_102):  # 103 colors
                palette = default_102
            else:
                rgb_rainbow = mpl.cm.rainbow(np.linspace(0,1,length))
                palette = [mpl.colors.rgb2hex(rgb_rainbow[i,:-1]) for i in range(length)]
        for i,x in enumerate(sorted(categories)):
            colors.append(palette[i])
        colors = dict(zip(categories, colors))
    del adata
    gc.collect()
    return jsonify({ann: colors, "labels": categories})


def get_paga3d_pos(adata):
    assert (
        adata.obsm["X_umap"].shape[1] >= 3
    ), """The embedding space should have at least three dimensions. 
        please set `n_component = 3` in `sc.tl.umap()`"""
    groups = adata.obs[adata.uns["paga"]["groups"]]
    connectivities_coarse = adata.uns["paga"]["connectivities"]
    paga3d_pos = np.zeros((connectivities_coarse.shape[0], 3))
    for i in range(connectivities_coarse.shape[0]):
        subset = (groups == groups.cat.categories[i]).values
        paga3d_pos[i] = np.median(adata.obsm["X_umap"][subset], axis=0)
    return paga3d_pos


@app.server.route("/columns", methods=["GET", "POST"])
def get_available_annotations():
    """
    http://127.0.0.1:8000/columns?db_name=1_scanpy_10xpbmc
    """
    db_name = request.args.get("db_name")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0]

    try:
        del adata
        gc.collect()
    except:
        pass

    adata = None
    if get_dataset_type_adata(db_name).lower() in ["scanpy", "velocity", "seurat", "paga"]:
        adata = sc.read(filename)
    else:
        adata = st.read(filename, file_format="pkl", workdir="./")
    
    annotations = [name for name in list(adata.obs.columns) if name not in ['branch_id', 'branch_id_alias']]
    del adata
    gc.collect()
    # Hack to remove two stream adata annotations that dont work in the annotation menu
    return jsonify(annotations)


def get_available_annotations_adata(adata):
    return adata.obs.columns


@app.server.route("/genes", methods=["GET", "POST"])
def get_genes():
    """
    http://127.0.0.1:8000/genes?db_name=1_scanpy_10xpbmc
    """
    db_name = request.args.get("db_name")
    try:
        del adata
        gc.collect()
    except:
        pass
    adata = None
    if get_dataset_type_adata(db_name).lower() == 'stream':
        adata = st.read(glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0], file_format="pkl", workdir="./")
    else:
        adata = sc.read(glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0])
    
    genes = adata.var_names
    if adata:
        del adata
        gc.collect()
    return jsonify(list(genes))


def get_genes_adata(adata):
    return adata.var_names


@app.server.route("/ts", methods=["GET", "POST"])
def get_ts():
    """
    velocity examples:
    http://127.0.0.1:8000/ts?db_name=3_velocity_pancrease&feature=clusters
    """

    db_name = request.args.get("db_name")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0]

    try:
        del adata
    except:
        pass

    if get_dataset_type_adata(db_name) == "velocity":
        adata = sc.read(filename)
        ts = [k.replace('absolute_velocity_umap_', '').replace('s', '')
              for k in adata.obsm.keys() if k.startswith('absolute')]
    del adata
    gc.collect()
    return jsonify(list(ts))

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
    return datasets_payload()


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
    img = qrcode.make(API + "/view/" + str(unique_id))
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
        else:
            file_id = value + '_scvr'
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
def update_output(unique_id, file_id):
    if not unique_id and not file_id:
        return dbc.Button("Please choose or upload dataset!", id='button', className="fly-button", color="link", disabled=True,n_clicks=0)
    if 'scvr' in file_id:
        return html.A(dbc.Button("Let's fly!", id='button',disabled=False,n_clicks=0, color="link", className="fly-button"),
                                    href="/view?dataset="+str(file_id) + "&fulldataset=false")
    else:     
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


### modifed from scanpy palettes https://github.com/theislab/scanpy/blob/master/scanpy/plotting/palettes.py

"""Color palettes in addition to matplotlib's palettes."""

from matplotlib import cm, colors

# Colorblindness adjusted vega_10
# See https://github.com/theislab/scanpy/issues/387
vega_10 = list(map(colors.to_hex, cm.tab10.colors))
vega_10_scanpy = vega_10.copy()
vega_10_scanpy[2] = '#279e68'  # green
vega_10_scanpy[4] = '#aa40fc'  # purple
vega_10_scanpy[8] = '#b5bd61'  # kakhi

# default matplotlib 2.0 palette
# see 'category20' on https://github.com/vega/vega/wiki/Scales#scale-range-literals
vega_20 = list(map(colors.to_hex, cm.tab20.colors))

# reorderd, some removed, some added
vega_20_scanpy = [
    *vega_20[0:14:2], *vega_20[16::2],  # dark without grey
    *vega_20[1:15:2], *vega_20[17::2],  # light without grey
    '#ad494a', '#8c6d31',  # manual additions
]
vega_20_scanpy[2] = vega_10_scanpy[2]
vega_20_scanpy[4] = vega_10_scanpy[4]
vega_20_scanpy[7] = vega_10_scanpy[8]  # kakhi shifted by missing grey
# TODO: also replace pale colors if necessary

default_20 = vega_20_scanpy

# https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
# update 1
# orig reference http://epub.wu.ac.at/1692/1/document.pdf
zeileis_28 = [
    "#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3",
    "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593",
    "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7",
    "#f3e1eb", "#f6c4e1", "#f79cd4",
    '#7f7f7f', "#c7c7c7", "#1CE6FF", "#336600",  # these last ones were added,
]

default_28 = zeileis_28

# from http://godsnotwheregodsnot.blogspot.de/2012/09/color-distribution-methodology.html
godsnot_102 = [
    # "#000000",  # remove the black, as often, we have black colored annotation
    "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72",
]

default_102 = godsnot_102

if __name__ == "__main__":
    # app.run_server(port=os.environ['PORT'])
    app.run_server(port=8050,debug=True,host='0.0.0.0',ssl_context='adhoc')

