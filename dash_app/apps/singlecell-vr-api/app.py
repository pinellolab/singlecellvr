import scanpy as sc
from scipy.sparse import isspmatrix
from glob import glob
import pathlib
import os
from flask import Flask, jsonify, request
import matplotlib as mpl
import networkx as nx
import numpy as np
#import stream as st
import pandas as pd


# Initialize app
APP_PATH = str(pathlib.Path(__file__).parent.resolve())
DATASET_DIRECTORY = os.path.join(APP_PATH, "app_datasets")
UPLOAD_DIRECTORY = os.path.join(APP_PATH, "app_uploaded_files")
QR_DIRECTORY = os.path.join(APP_PATH, "assets")
server = Flask(__name__)
API = "http://localhost:8000/"

if not os.path.exists(DATASET_DIRECTORY):
    os.makedirs(DATASET_DIRECTORY)


@server.route("/databases", methods=["GET", "POST"])
def get_databases():
    adata_files = glob(os.path.join(DATASET_DIRECTORY, "*"))
    adata_list = list()
    for i in adata_files:
        label = os.path.basename(i).split(".")[0]
        label_list = label.split("_")
        # ID: [LABEL, PATH]
        adata_list.append(
            {
                "value": label_list[0],
                "label": label,
                "path": i,
                "type": label_list[1],
            }
        )
    return jsonify(adata_list)


@server.route("/data_type", methods=["GET", "POST"])
def get_dataset_type():
    """
    http://127.0.0.1:8000/data_type?db_name=1_scanpy_10xpbmc
    """
    db_name = request.args.get("db_name")
    return jsonify({"type": db_name.split("_")[1]})


def get_dataset_type_adata(db_name):
    return db_name.split("_")[1]


@server.route("/coordinates", methods=["GET", "POST"])
def get_cooridinates():
    """
    http://127.0.0.1:8000/coordinates?db_name=1_scanpy_10xpbmc&embed=umap
    http://127.0.0.1:8000/coordinates?db_name=3_velocity_pancrease&embed=umap
    http://127.0.0.1:8000/coordinates?db_name=4_seurat_10xpbmc&embed=umap
    http://127.0.0.1:8000/coordinates?db_name=5_stream_nestorowa16&embed=umap
    """
    db_name = request.args.get("db_name")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0]

    if get_dataset_type_adata(db_name) in ["scanpy", "velocity", "seurat"]:
        adata = sc.read(filename)
        embed = request.args.get("embed")
    else:
        adata = st.read(filename, file_format="pkl", workdir="./")

    list_cells = []
    for i in range(adata.shape[0]):
        dict_coord_cells = dict()
        dict_coord_cells["cell_id"] = adata.obs_names[i]
        if get_dataset_type_adata(db_name) in ["scanpy", "velocity"]:
            dict_coord_cells["x"] = str(adata.obsm[f"X_{embed}"][i, 0])
            dict_coord_cells["y"] = str(adata.obsm[f"X_{embed}"][i, 1])
            dict_coord_cells["z"] = str(adata.obsm[f"X_{embed}"][i, 2])
        elif get_dataset_type_adata(db_name) == "seurat":
            dict_coord_cells["x"] = str(adata.obsm[f"{embed}_cell_embeddings"][i, 0])
            dict_coord_cells["y"] = str(adata.obsm[f"{embed}_cell_embeddings"][i, 1])
            dict_coord_cells["z"] = str(adata.obsm[f"{embed}_cell_embeddings"][i, 2])
        elif get_dataset_type_adata(db_name) == "stream":
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

            dict_coord_cells["x"] = str(adata.obsm["X_dr"][i, 0])
            dict_coord_cells["y"] = str(adata.obsm["X_dr"][i, 1])
            dict_coord_cells["z"] = str(adata.obsm["X_dr"][i, 2])
            return jsonify(
                {"nodes": dict_nodes, "edges": list_edges, "graph": list_curves}
            )
        else:
            raise TypeError("not supported format")
        list_cells.append(dict_coord_cells)
    return jsonify(list_cells)


@server.route("/features", methods=["GET", "POST"])
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
    if feature == "velocity":
        embed = request.args.get("embed")

    if get_dataset_type_adata(database) in ["scanpy", "velocity", "seurat"]:
        adata = sc.read(filename)
    else:
        adata = st.read(filename, file_format="pkl", workdir="./")

    list_metadata = []
    if feature in get_available_annotations_adata(adata):  # cluster columns
        dict_colors = {
            feature: dict(
                zip(adata.obs[feature].cat.categories, adata.uns[f"{feature}_colors"])
            )
        }
        for i in range(adata.shape[0]):
            dict_metadata = dict()
            dict_metadata["cell_id"] = adata.obs_names[i]
            dict_metadata["clusters"] = adata.obs[feature].tolist()[i]
            dict_metadata["clusters_color"] = dict_colors[feature][
                dict_metadata["clusters"]
            ]
            list_metadata.append(dict_metadata)
    elif feature in ["expression", "rna"] or feature in get_available_annotations_adata(
        adata
    ):  # pseudotime or latent_time columns
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

            dict_coord_cells["x0"] = str(adata.obsm[f"X_{embed}"][i, 0])
            dict_coord_cells["y0"] = str(adata.obsm[f"X_{embed}"][i, 1])
            dict_coord_cells["z0"] = str(adata.obsm[f"X_{embed}"][i, 2])

            if time == "None":
                dict_coord_cells["x1"] = str(adata.obsm[f"velocity_{embed}"][i, 0])
                dict_coord_cells["y1"] = str(adata.obsm[f"velocity_{embed}"][i, 1])
                dict_coord_cells["z1"] = str(adata.obsm[f"velocity_{embed}"][i, 2])
            elif time in list(map(str, [0.1, 1, 5, 10, 20, 30, 50, 100])):
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
            #if isinstance(adata.obs_names[i], bytes):
            #    dict_coord_cells["cell_id"] = adata.obs_names[i].decode("utf-8")
            #else:
            #    dict_coord_cells["cell_id"] = adata.obs_names[i]

            if time == "None":
                dict_coord_cells["x0"] = str(adata.uns[f"X_grid"][i, 0])
                dict_coord_cells["y0"] = str(adata.uns[f"X_grid"][i, 1])
                dict_coord_cells["z0"] = str(adata.uns[f"X_grid"][i, 2])
                dict_coord_cells["x1"] = str(adata.uns[f"V_grid"][i, 0])
                dict_coord_cells["y1"] = str(adata.uns[f"V_grid"][i, 1])
                dict_coord_cells["z1"] = str(adata.uns[f"V_grid"][i, 2])
            elif time in list(map(str, [0.01, 0.1, 1, 5, 10, 20, 50, 80, 100])):
                dict_coord_cells["x0"] = str(adata.uns[f"X_grid_{time}"][i, 0])
                dict_coord_cells["y0"] = str(adata.uns[f"X_grid_{time}"][i, 1])
                dict_coord_cells["z0"] = str(adata.uns[f"X_grid_{time}"][i, 2])
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
    return jsonify(list_metadata)


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


@server.route("/columns", methods=["GET", "POST"])
def get_available_annotations():
    """
    http://127.0.0.1:8000/columns?db_name=1_scanpy_10xpbmc
    """
    db_name = request.args.get("db_name")
    adata = sc.read(glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0])
    return jsonify(list(adata.obs.columns))


def get_available_annotations_adata(adata):
    return adata.obs.columns


@server.route("/genes", methods=["GET", "POST"])
def get_genes():
    """
    http://127.0.0.1:8000/genes?db_name=1_scanpy_10xpbmc
    """
    db_name = request.args.get("db_name")
    adata = sc.read(glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0])
    return jsonify(list(adata.var_names))


def get_genes_adata(adata):
    return adata.var_names


@server.route("/ts", methods=["GET", "POST"])
def get_ts():
    """
    velocity examples:
    http://127.0.0.1:8000/ts?db_name=3_velocity_pancrease
    """

    db_name = request.args.get("db_name")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0]

    if get_dataset_type_adata(db_name) == "velocity":
        adata = sc.read(filename)
        ts = [k.replace('absolute_velocity_umap_', '').replace('s', '')
              for k in adata.obsm.keys() if k.startswith('absolute')]
    return jsonify(list(ts))
