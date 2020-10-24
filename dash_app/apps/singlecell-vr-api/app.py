import scanpy as sc
from scipy.sparse import isspmatrix

from glob import glob
import pathlib
import os
from flask import Flask, jsonify, request
import matplotlib as mpl

# Initialize app
APP_PATH = str(pathlib.Path(__file__).parent.resolve())
DATASET_DIRECTORY = os.path.join(APP_PATH, "app_datasets")
UPLOAD_DIRECTORY = os.path.join(APP_PATH, "app_uploaded_files")
QR_DIRECTORY = os.path.join(APP_PATH, "assets")
server = Flask(__name__)


@server.route('/databases', methods=['GET', 'POST'])
def get_databases():
    adata_files = glob(os.path.join(DATASET_DIRECTORY, f'*'))
    adata_list = list()
    for i in adata_files:
        label = os.path.basename(i).split('.')[0]
        label_list = label.split('_')
        label = ' '.join(label_list[1:])
        # ID: [LABEL, PATH]
        adata_list.append({'value': label_list[0], 'label': label.upper(), 'path': i})
    return jsonify(adata_list)

@server.route('/coordinates', methods=['GET', 'POST'])
def get_cooridinates():
    database = request.args.get('database')
    adata = sc.read(glob(os.path.join(DATASET_DIRECTORY, f'{database}_*'))[0])
    list_cells = []
    for i in range(adata.shape[0]):
        dict_coord_cells = dict()
        dict_coord_cells['cell_id'] = adata.obs_names[i]
        dict_coord_cells['x'] = str(adata.obsm['X_umap'][i, 0])
        dict_coord_cells['y'] = str(adata.obsm['X_umap'][i, 1])
        dict_coord_cells['z'] = str(adata.obsm['X_umap'][i, 2])
        list_cells.append(dict_coord_cells)
    return jsonify(list_cells)

@server.route('/features', methods=['GET', 'POST'])
def get_features():
    database = request.args.get('database')
    feature = request.args.get('feature')

    adata = sc.read(glob(os.path.join(DATASET_DIRECTORY, f'{database}_*'))[0])
    list_metadata = []
    if feature == 'clusters':
        dict_colors = {'clusters': dict(zip(adata.obs['clusters'].cat.categories,
                                            adata.uns['clusters_colors']))}
        for i in range(adata.shape[0]):
            dict_metadata = dict()
            dict_metadata['cell_id'] = adata.obs_names[i]
            dict_metadata['clusters'] = adata.obs['clusters'].tolist()[i]
            dict_metadata['clusters_color'] = dict_colors['clusters'][dict_metadata['clusters']]
            list_metadata.append(dict_metadata)
    elif feature == 'expression':
        gene = request.args.get('gene')
        if gene not in adata.var_names:
            return jsonify({})
        else:
            expr = adata[:, gene].X.toarray()[:, 0] if isspmatrix(adata.X) else adata[:, gene].X
            cm = mpl.cm.get_cmap('viridis', 512)
            norm = mpl.colors.Normalize(vmin=0, vmax=max(expr), clip=True)
            list_metadata = []
            for i, x in enumerate(adata.obs_names):
                dict_genes = dict()
                dict_genes['cell_id'] = x
                dict_genes['color'] = mpl.colors.to_hex(cm(norm(expr[i])))
                list_metadata.append(dict_genes)
    return jsonify(list_metadata)
