import numpy as np
import pandas as pd
import os
import json
import shutil
import networkx as nx
import matplotlib as mpl
from scipy.sparse import isspmatrix
from pandas.api.types import is_string_dtype,is_numeric_dtype
# import scvelo as scv

from . import palettes

def get_colors(adata,ann):
    df_cell_colors = pd.DataFrame(index=adata.obs.index)
    if(is_numeric_dtype(adata.obs[ann])):
        cm = mpl.cm.get_cmap('viridis',512)
        norm = mpl.colors.Normalize(vmin=0, vmax=max(adata.obs[ann]),clip=True)
        df_cell_colors[ann+'_color'] = [mpl.colors.to_hex(cm(norm(x))) for x in adata.obs[ann]]
    else:
        df_cell_colors[ann+'_color'] = ''

        adata.obs[ann] = adata.obs[ann].astype('category')
        categories = adata.obs[ann].cat.categories
        length = len(categories)
        # check if default matplotlib palette has enough colors
        # mpl.style.use('default')
        if len(mpl.rcParams['axes.prop_cycle'].by_key()['color']) >= length:
            cc = mpl.rcParams['axes.prop_cycle']()
            palette = [next(cc)['color'] for _ in range(length)]
        else:
            if length <= 20:
                palette = palettes.default_20
            elif length <= 28:
                palette = palettes.default_28
            elif length <= len(palettes.default_102):  # 103 colors
                palette = palettes.default_102
            else:
                rgb_rainbow = mpl.cm.rainbow(np.linspace(0,1,length))
                palette = [mpl.colors.rgb2hex(rgb_rainbow[i,:-1]) for i in range(length)]
        for i,x in enumerate(categories):
            id_cells = np.where(adata.obs[ann]==x)[0]
            df_cell_colors.loc[df_cell_colors.index[id_cells],ann+'_color'] = palette[i]
    return(df_cell_colors[ann+'_color'].tolist())

def output_scanpy_cells(adata,ann_list,reportdir='./scanpy_report',gene_list=None):
    assert (adata.obsm['X_umap'].shape[1]>=3),\
    '''The embedding space should have at least three dimensions. 
    please set `n_component = 3` in `sc.tl.umap()`'''
    ###remove duplicate keys
    ann_list = list(dict.fromkeys(ann_list))
    ### make sure all labels exist
    for ann in ann_list:
        if ann not in adata.obs.columns:
            raise ValueError('could not find %s in %s'  % (ann,adata.obs.columns))
    if(gene_list is not None):
        ###remove duplicate keys
        gene_list = list(dict.fromkeys(gene_list))
        for gene in gene_list:
            if(gene not in adata.raw.var_names):
                raise ValueError('could not find %s in `adata.var_names`'  % (gene))
    try:
        if(not os.path.exists(reportdir)):
                os.makedirs(reportdir)
        ## output coordinates of cells
        list_cells = []
        for i in range(adata.shape[0]):
            dict_coord_cells = dict()
            dict_coord_cells['cell_id'] = adata.obs_names[i]
            dict_coord_cells['x'] = str(adata.obsm['X_umap'][i,0])
            dict_coord_cells['y'] = str(adata.obsm['X_umap'][i,1])
            dict_coord_cells['z'] = str(adata.obsm['X_umap'][i,2])
            list_cells.append(dict_coord_cells)
        with open(os.path.join(reportdir,'scatter.json'), 'w') as f:
            json.dump(list_cells, f)

        ## output metadata file of cells
        list_metadata = []
        dict_colors = get_paga_colors(adata,ann_list)
        for i in range(adata.shape[0]):
            dict_metadata = dict()
            dict_metadata['cell_id'] = adata.obs_names[i]
            for ann in ann_list:
                dict_metadata[ann] = adata.obs[ann].tolist()[i]
                dict_metadata[ann+'_color'] = dict_colors[ann][i]
            list_metadata.append(dict_metadata)
        with open(os.path.join(reportdir,'metadata.json'), 'w') as f:
            json.dump(list_metadata, f)

        ## output gene expression of cells
        if(gene_list is not None):
            df_genes = pd.DataFrame(adata.raw.X.toarray() if isspmatrix(adata.raw.X) else adata.raw.X,
                                    index=adata.raw.obs_names,
                                    columns=adata.raw.var_names)
            cm = mpl.cm.get_cmap('viridis',512)
            for g in gene_list:
                list_genes = []
                norm = mpl.colors.Normalize(vmin=0, vmax=max(df_genes[g]),clip=True)
                for x in adata.obs_names:
                    dict_genes = dict()
                    dict_genes['cell_id'] = x
                    dict_genes['color'] = mpl.colors.to_hex(cm(norm(df_genes.loc[x,g])))
                    list_genes.append(dict_genes)
                with open(os.path.join(reportdir,'gene_'+g+'.json'), 'w') as f:
                    json.dump(list_genes, f)
    except:
        print("Output cells: failed!")
        raise
    else:
        print("Output cells: finished!")

def get_paga_colors(adata,ann_list):
    dict_colors = dict()
    for ann in ann_list:
        if(ann+'_colors' in adata.uns_keys()):
            ### if ann can be found in adata, then the same colors will be used
            df_cell_colors = pd.DataFrame(index=adata.obs.index)
            df_cell_colors[ann+'_color'] = ''
            adata.obs[ann] = adata.obs[ann].astype('category')
            categories = adata.obs[ann].cat.categories
            for i,x in enumerate(categories):
                id_cells = np.where(adata.obs[ann]==x)[0]
                df_cell_colors.loc[df_cell_colors.index[id_cells],ann+'_color'] = adata.uns[ann+'_colors'][i]
            dict_colors[ann] = df_cell_colors[ann+'_color'].tolist()
        else:
            dict_colors[ann] = get_colors(adata,ann)
    return(dict_colors)

def get_paga3d_pos(adata):
    assert (adata.obsm['X_umap'].shape[1]>=3),\
    '''The embedding space should have at least three dimensions. 
    please set `n_component = 3` in `sc.tl.umap()`'''
    groups = adata.obs[adata.uns['paga']['groups']]
    connectivities_coarse = adata.uns['paga']['connectivities']
    paga3d_pos = np.zeros((connectivities_coarse.shape[0], 3))
    for i in range(connectivities_coarse.shape[0]):
        subset = (groups == groups.cat.categories[i]).values
        paga3d_pos[i] = np.median(adata.obsm['X_umap'][subset],axis=0)
    return paga3d_pos

def output_paga_graph(adata,node_name = None,reportdir='./paga_report'):
    try:
        if(not os.path.exists(reportdir)):
                os.makedirs(reportdir)
        G = nx.from_numpy_matrix(adata.uns['paga']['connectivities'].toarray())
        adata.uns['paga']['pos'] = get_paga3d_pos(adata)
        ## output coordinates of paga graph
        list_lines = []
        for edge_i in G.edges():
            dict_coord_lines = dict()
            dict_coord_lines['branch_id'] = [[str(edge_i[0]),str(edge_i[1])]]
            dict_coord_lines['xyz'] = [{'x':pos[0],'y':pos[1],'z':pos[2]} for pos in adata.uns['paga']['pos'][[edge_i[0],edge_i[1]],:]]
            list_lines.append(dict_coord_lines)
        with open(os.path.join(reportdir,'graph_paths.json'), 'w') as f:
            json.dump(list_lines, f)

        ## output topology of paga graph
        dict_nodes = dict()
        list_edges = []
        if(node_name is None):
            dict_nodename = {i: adata.obs[adata.uns['paga']['groups']].cat.categories[i] for i in G.nodes()}
        else:
            if(isinstance(node_name,dict)):
                dict_nodename = node_name
            else:
                raise TypeError('node_name must be a dict')
        for node_i in G.nodes():
            dict_nodes_i = dict()
            dict_nodes_i['node_name'] = dict_nodename[node_i]
            dict_nodes_i['xyz'] = {'x':adata.uns['paga']['pos'][:,0][node_i],
                                   'y':adata.uns['paga']['pos'][:,1][node_i],
                                   'z':adata.uns['paga']['pos'][:,2][node_i]}
            dict_nodes[node_i] = dict_nodes_i
        for edge_i in G.edges():
            dict_edges = dict()
            dict_edges['nodes'] = [str(edge_i[0]),str(edge_i[1])]
            dict_edges['weight'] = adata.uns['paga']['connectivities'][edge_i[0],edge_i[1]]
            list_edges.append(dict_edges)
        with open(os.path.join(reportdir,'graph_nodes.json'), 'w') as f:
            json.dump(dict_nodes, f)
        with open(os.path.join(reportdir,'graph_edges.json'), 'w') as f:
            json.dump(list_edges, f)
    except:
        print("Output graph: failed!")
        raise
    else:
        print("Output graph: finished!")

def output_paga_cells(adata,ann_list,reportdir='./paga_report',gene_list=None):
    output_scanpy_cells(adata,ann_list,reportdir=reportdir,gene_list=gene_list)


def output_seurat_cells(adata,ann_list,reportdir='./seurat_report',gene_list=None):
    assert (adata.obsm['umap_cell_embeddings'].shape[1]>=3),\
    '''The embedding space should have at least three dimensions. 
    please set `n.component = 3` in `RunUMAP()`'''
    ###remove duplicate keys
    ann_list = list(dict.fromkeys(ann_list))
    ### make sure all labels exist
    for ann in ann_list:
        ann = ann.strip(" ").strip("\n")
        if ann not in adata.obs.columns:
            raise ValueError('could not find %s in %s'  % (ann,adata.obs.columns))
    if(gene_list is not None):
        ###remove duplicate keys
        gene_list = list(dict.fromkeys(gene_list))
        for gene in gene_list:
            gene = gene.strip(" ").strip("\n")
            if(gene not in adata.var_names):
                raise ValueError('could not find %s in `adata.var_names`'  % (gene))
    try:
        if(not os.path.exists(reportdir)):
                os.makedirs(reportdir)
        ## output coordinates of cells
        list_cells = []
        for i in range(adata.shape[0]):
            dict_coord_cells = dict()
            dict_coord_cells['cell_id'] = adata.obs_names[i]
            dict_coord_cells['x'] = str(adata.obsm['umap_cell_embeddings'][i,0])
            dict_coord_cells['y'] = str(adata.obsm['umap_cell_embeddings'][i,1])
            dict_coord_cells['z'] = str(adata.obsm['umap_cell_embeddings'][i,2])
            list_cells.append(dict_coord_cells)
        with open(os.path.join(reportdir,'scatter.json'), 'w') as f:
            json.dump(list_cells, f)

        ## output metadata file of cells
        list_metadata = []
        dict_colors = dict()
        for ann in ann_list:
            dict_colors[ann] = get_colors(adata,ann)
        for i in range(adata.shape[0]):
            dict_metadata = dict()
            dict_metadata['cell_id'] = adata.obs_names[i]
            for ann in ann_list:
                dict_metadata[ann] = adata.obs[ann].tolist()[i]
                dict_metadata[ann+'_color'] = dict_colors[ann][i]
            list_metadata.append(dict_metadata)
        with open(os.path.join(reportdir,'metadata.json'), 'w') as f:
            json.dump(list_metadata, f)

        ## output gene expression of cells
        if(gene_list is not None):
            df_genes = pd.DataFrame(adata.layers['norm_data'].toarray() if isspmatrix(adata.layers['norm_data']) else adata.layers['norm_data'],
                                    index=adata.obs_names,
                                    columns=adata.var_names)
            cm = mpl.cm.get_cmap('viridis',512)
            for g in gene_list:
                list_genes = []
                norm = mpl.colors.Normalize(vmin=0, vmax=max(df_genes[g]),clip=True)
                for x in adata.obs_names:
                    dict_genes = dict()
                    dict_genes['cell_id'] = x
                    dict_genes['color'] = mpl.colors.to_hex(cm(norm(df_genes.loc[x,g])))
                    list_genes.append(dict_genes)
                with open(os.path.join(reportdir,'gene_'+g+'.json'), 'w') as f:
                    json.dump(list_genes, f)
    except:
        print("Output cells: failed!")
        raise
    else:
        print("Output cells: finished!")


def output_velocity_cells(adata, ann_field, gene_list=None,
                          reportdir='./velocity_report'):

    if gene_list is not None:
        ###remove duplicate keys
        gene_list = list(dict.fromkeys(gene_list))
        for gene in gene_list:
            if gene not in adata.raw.var_names:
                raise ValueError('could not find %s in `adata.var_names`'  % (gene))
    try:
        if not os.path.exists(reportdir):
            os.makedirs(reportdir)

        list_cells = []
        for i in range(adata.shape[0]):
            dict_coord_cells = dict()
            dict_coord_cells['cell_id'] = adata.obs_names[i]
            dict_coord_cells['x0'] = str(adata.obsm['X_umap'][i,0])
            dict_coord_cells['y0'] = str(adata.obsm['X_umap'][i,1])
            dict_coord_cells['z0'] = str(adata.obsm['X_umap'][i,2])

            dict_coord_cells['x1'] = str(adata.obsm['velocity_umap'][i,0])
            dict_coord_cells['y1'] = str(adata.obsm['velocity_umap'][i,1])
            dict_coord_cells['z1'] = str(adata.obsm['velocity_umap'][i,2])
            list_cells.append(dict_coord_cells)
        with open(os.path.join(reportdir, 'scatter.json'), 'w') as f:
            json.dump(list_cells, f)

        list_metadata = []
        print(adata.uns)

        dict_colors = {ann_field: dict(zip(adata.obs[ann_field].cat.categories,
                                           adata.uns[f'{ann_field}_colors']))}
        print(dict_colors)
        for i in range(adata.shape[0]):
            dict_metadata = dict()
            dict_metadata['cell_id'] = adata.obs_names[i]
            dict_metadata[ann_field] = adata.obs[ann_field].tolist()[i]
            dict_metadata[ann_field+'_color'] = dict_colors[ann_field][dict_metadata[ann_field]]
            list_metadata.append(dict_metadata)
        with open(os.path.join(reportdir,'metadata.json'), 'w') as f:
            json.dump(list_metadata, f)

        if gene_list is not None:
            df_genes = pd.DataFrame(adata.raw.X.toarray() if isspmatrix(adata.raw.X) else adata.raw.X,
                                    index=adata.raw.obs_names,
                                    columns=adata.raw.var_names)
            cm = mpl.cm.get_cmap('viridis',512)
            for g in gene_list:
                list_genes = []
                norm = mpl.colors.Normalize(vmin=0, vmax=max(df_genes[g]),clip=True)
                for x in adata.obs_names:
                    dict_genes = dict()
                    dict_genes['cell_id'] = x
                    dict_genes['color'] = mpl.colors.to_hex(cm(norm(df_genes.loc[x,g])))
                    list_genes.append(dict_genes)
                with open(os.path.join(reportdir,'gene_'+g+'.json'), 'w') as f:
                    json.dump(list_genes, f)
    except:
        print("Output velocity failed")