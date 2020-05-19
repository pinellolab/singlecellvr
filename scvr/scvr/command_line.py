#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Authors: Huidong Chen
# Contact information: huidong.chen@mgh.harvard.edu

import warnings
warnings.filterwarnings('ignore')

__tool_name__='scvr'

import pandas as pd
import anndata as ad
import json
import argparse
import os
import shutil
import scvr

print('- Single cell VR preprocessing -',flush=True)
print('Version %s\n' % scvr.__version__,flush=True)


def main():
    parser = argparse.ArgumentParser(description='%s Parameters' % __tool_name__ ,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--filename", dest="filename",default = None,required=True,
                        help="Analysis result file name", metavar="FILE")
    parser.add_argument("-t", "--toolname",dest="toolname", default=None,required=True,
                        type = str.lower,choices=['scanpy','paga','seurat','stream'],
                        help="Tool used to generate the analysis result.")
    parser.add_argument("-a","--annotations",dest="annotations", default=None,required=True,
                        help="Annotation file name. It contains the cell annotation key(s) to visualize in one column.")
    parser.add_argument("-g","--genes",dest="genes", default=None,
                        help="Gene list file name. It contains the genes to visualize in one column.")
    parser.add_argument("-o","--output",dest="output", default='vr_report',
                        help="Output folder name")

    args = parser.parse_args()
    filename = args.filename
    toolname = args.toolname
    genes = args.genes
    output = args.output #work directory
    annotations = args.annotations

    if(annotations is None):
        raise Exception("Annotation file must be specified when %s is chosen." % (toolname))
    try:
        ann_list = pd.read_csv(annotations,sep='\t',header=None,index_col=None).iloc[:,0].tolist()
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise
    except:
        print('Failed to load in annotation file.')
        raise
    else:
        ann_list = list(set(ann_list))

    if(genes is not None):
        try:
            gene_list = pd.read_csv(genes,sep='\t',header=None,index_col=None).iloc[:,0].tolist()
        except FileNotFoundError as fnf_error:
            print(fnf_error)
            raise
        except:
            print('Failed to load in gene list.')
            raise
        else:
            gene_list = list(set(gene_list))
    else:
        gene_list = None

    print("Converting '%s' analysis result ..." % toolname)
    
    if(toolname in ['scanpy','paga','seurat']):
        if(toolname=='scanpy'):
            assert (filename.lower().endswith(('.h5ad'))), "For PAGA only .h5ad file is supported."
            print('reading in h5ad file ...')
            adata = ad.read_h5ad(filename)
            scvr.output_scanpy_cells(adata,ann_list,gene_list=gene_list,reportdir=output)
        if(toolname=='paga'):
            assert (filename.lower().endswith(('.h5ad'))), "For PAGA only .h5ad file is supported."
            print('reading in h5ad file ...')
            adata = ad.read_h5ad(filename)
            scvr.output_paga_graph(adata,reportdir=output)
            scvr.output_paga_cells(adata,ann_list,gene_list=gene_list,reportdir=output)
        if(toolname=='seurat'):
            assert (filename.lower().endswith(('.loom'))), "For Seurat only .loom file is supported."
            print('reading in loom file ...')
            adata = ad.read_loom(filename)
            scvr.output_seurat_cells(adata,ann_list,gene_list=gene_list,reportdir=output)
        with open(os.path.join(output,'index.json'), 'w') as f:
                json.dump({ "tool": toolname }, f)
        shutil.make_archive(base_name=output, format='zip',root_dir=output)
        shutil.rmtree(output)
    if(toolname=='stream'):
        try:
            import stream as st
        except ImportError:
            raise ImportError(
                'Please install STREAM >=0.5: `conda install -c bioconda stream`.'
            )
        assert (filename.lower().endswith(('.pkl'))), "For STREAM only .pkl file is supported."
        print('reading in pkl file ...')
        adata = st.read(filename,file_format='pkl',workdir='./')
        st.save_vr_report(adata,ann_list=ann_list,gene_list=gene_list,file_name=output)

if __name__ == "__main__":
    main()
