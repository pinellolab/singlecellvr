#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Authors: Huidong Chen
# Contact information: huidong.chen@mgh.harvard.edu

import warnings
warnings.filterwarnings('ignore')

__tool_name__='scvr_prep'

import pandas as pd
import argparse
import os
import shutil
import scvr_prep

print('- Single cell VR preprocessing -',flush=True)
print('Version %s\n' % scvr_prep.__version__,flush=True)


def main():
    parser = argparse.ArgumentParser(description='%s Parameters' % __tool_name__ ,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--filename", dest="filename",default = None,required=True,
                        help="analysis result file name", metavar="FILE")
    parser.add_argument("-t", "--toolname",dest="toolname", default=None,required=True,
                        help="tool name used to generate the analysis result")
    parser.add_argument("-g","--genes",dest="genes", default=None,
                        help="genes to visualize. it can either be filename which contains all the genes in one column or a set of gene names separated by comma")
    parser.add_argument("-l","--label",dest="label", default=None,
                        help="label used to color cells")
    parser.add_argument("-o","--output",dest="output", default='vr_report',
                        help="output folder")

    args = parser.parse_args()

    filename = args.filename
    toolname = args.toolname
    genes = args.genes
    output = args.output #work directory
    label = args.label
    if(os.path.exists(genes)):
        gene_list = pd.read_csv(genes,sep='\t',header=None,index_col=None).iloc[:,0].tolist()
        gene_list = list(set(gene_list))
    else:
        gene_list = genes.split(',')
    if(toolname=='paga'):
        import scanpy as sc
        adata = sc.read_h5ad(filename)
        adata.uns['paga']['pos'] = scvr_prep.get_paga3d_pos(adata)
        scvr_prep.output_paga_graph(adata,reportdir=output)
        scvr_prep.output_paga_cells(adata,label=label,genes=gene_list,reportdir=output)
        shutil.make_archive(base_name=output, format='zip',root_dir=output)
        shutil.rmtree(output)
    if(toolname=='stream'):
        import stream as st
        adata = st.read(filename,file_format='pkl',workdir='./')
        st.save_vr_report(adata,genes=gene_list,file_name=output)

if __name__ == "__main__":
    main()