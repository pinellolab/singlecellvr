# singlecellvr

<img src="images/SCVR_logo.png" alt="http://www.singlecellvr.com" width="400" height="160">

<img src="images/scvr.jpeg" alt="http://www.singlecellvr.com" width="400" height="40">

Single cell visualization using Virtual Reality (VR)  

http://www.singlecellvr.com/

SingleCellVR can be used with our preprocessed datasets found at the link above or by following the steps below to process your own dataset.

## Running singlecellVR locally

singlecellVR is available to be run locally with Docker. The following steps demonstrate how to use singlecellVR with Docker.

### Step 1:
Create a directory and move your datasets into it. The datasets are your h5ad, loom, pkl, or scvr converted zip files containing your single-cell data. The datasets referenced in the singlecellVR manuscript are available for download here: https://www.dropbox.com/sh/4zoaost27ky91ob/AADMJXKdhgBpuO9qJfBgJ_V6a?dl=1.

### Step 2:
Run the docker image. In order to serve the singlecellVR webapp locally with Docker, you will need two free ports on your system. The docker run command requires
the following arguments:
  - {path}, the absolute path to the directory you created in step 1, e.g. "/Users/david/singlecell/"
  - {backend-port}, a free port on your system over which to serve the singlecellVR backend, e.g. 8000
  - {frontend-port}, a free port on your system over which to serve the singlecellVR frontend, e.g. 8080
  - {host-address}, the host address at which to serve singlecellVR, e.g. 0.0.0.0

The docker image may be run as follows 
`docker run --rm -v {path}:/app/dash_app/apps/dash-singlecell-vr/app_datasets/ -p {host-address}:{frontend-port}:{frontend-port} -p {host-address}:{backend-port}:{backend-port} pinellolab/singlecellvr {host-address} {frontend-port} {backend-port}`

The -v flag allows you to make your local datasets visible to the docker container. 

For example:
```
mkdir singlecell_data
cd singlecell_data
wget https://www.dropbox.com/sh/4zoaost27ky91ob/AADMJXKdhgBpuO9qJfBgJ_V6a?dl=1 -O singlecellvr_data.zip

# if you are on a Window machine extract the file manually or use this command line version:  wget http://stahlworks.com/dev/unzip.exe -O unzip.exe
./unzip singlecellvr_data.zip

# on windows you need to type instead:.\unzip singlecellvr_data.zip

rm singlecellvr_data.zip
docker run --rm -v ${PWD}:/app/dash_app/apps/dash-singlecell-vr/app_datasets/ -p 0.0.0.0:8080:8080 -p 0.0.0.0:8000:8000 pinellolab/singlecellvr https://0.0.0.0 8000 8080
```

Note: To enable VR capabilities the urls must be served over https. 

Now open the browser and open https://localhost:8080 to run a local version of singlecellvr on your machine

## SingleCellVR Preprocess:  

Prepare your data for the visualization on Single Cell VR website <https://singlecellvr.com/>

Installation
------------
Install and update using pip:  
`pip install scvr`

Usage
-----
`$ scvr --help`

```
usage: scvr [-h] -f FILE -t {scanpy,paga,seurat,stream} -a ANNOTATIONS [-g GENES] [-o OUTPUT] [--layer]

scvr Parameters

required arguments:
  -f FILE, --filename FILE
                        Analysis result file name (default: None)
  -t {scanpy,paga,seurat,stream}, --toolname {scanpy,paga,seurat,stream}
                        Tool used to generate the analysis result (default: None)
  -a ANNOTATIONS, --annotations ANNOTATIONS
                        Annotation file name. It contains the cell annotation key(s) 
                        to visualize in one column (default: None)
                        
optional arguments:

  -g GENES, --genes GENES
                        Gene list file name. It contains the genes 
                        to visualize in one column (default: None)
  -o OUTPUT, --output OUTPUT
                        Output folder name (default: scvr_report)
  --layer LAYER         The name of layer in Anndata object for gene
                        expression (default: norm_data)
  -h, --help            show this help message and exit
```


Examples:
---------
### Scanpy:  

To get single cell VR report for Scanpy :  
```bash
scvr -f ./scanpy_result/scanpy_10xpbmc.h5ad -t scanpy -a annotations.txt -g genes.txt -o scanpy_report
```

* Input files can be found [here](https://www.dropbox.com/sh/m6u9y38mi5qgf3o/AACe6cgiywaxM7ARtw54sg1Ha?dl=0) 
* To generate the `scanpy_10xpbmc.h5ad`, check out [Scanpy analysis](https://nbviewer.jupyter.org/github/pinellolab/singlecellvr/blob/master/examples/scanpy_10xpbmc.ipynb?flush_cache=true). *(Make sure set `n_components=3` in `sc.tl.umap(adata,n_components=3)`)*


### PAGA:  

To get single cell VR report for PAGA :  
```bash
scvr -f ./paga_result/paga3d_paul15.h5ad -t paga -a annotations.txt -g genes.txt -o paga_report
```

* Input files can be found [here](https://www.dropbox.com/sh/03zpxs9zv7yusi1/AADKVSU8Il1JcjA7lfHjmRpSa?dl=0) 
* To generate the `paga3d_paul15.h5ad`, check out [PAGA analysis](https://nbviewer.jupyter.org/github/pinellolab/singlecellvr/blob/master/examples/paga_paul15.ipynb?flush_cache=true). *(Make sure set `n_components=3` in `sc.tl.umap(adata,n_components=3)`)*

### Seurat:  
To get single cell VR report for Seurat scRNA-seq analysis:  
```bash
scvr -f ./seurat_result/seurat3d_10xpbmc.loom -t seurat -a annotations.txt -g genes.txt -o seurat_report
```
* Input files can be found [here](https://www.dropbox.com/sh/tpk4qfm5qsjpffn/AADmKmyDx7rhzKBOpIlAgMEUa?dl=0) 
* To generate the `seurat3d_10xpbmc.loom`, check out [Seurat analysis](https://nbviewer.jupyter.org/github/pinellolab/singlecellvr/blob/master/examples/seurat_10xpbmc.ipynb?flush_cache=true). *(Make sure set `n.components = 3` in `pbmc <- RunUMAP(pbmc, dims = 1:10, n.components = 3)`)*. NOTE the notebook is designed for Seurat 3.x based analysis, if you have installed Seurat 4.0.x, please use the [Seurat4 analysis](https://github.com/qinqian/singlecellvr/blob/master/examples/seurat4_pbmc3k.ipynb). 

To get single cell VR report for Seurat scRNA-seq integration analysis:
```bash
scvr -f ./seurat_integration_result/seurat3d_integration.loom -t seurat -a annotations.txt -g genes.txt -o seurat_integration_report --layer scale_data
```  
* To generate the `seurat3d_10xpbmc.loom`, check out [Seurat integration analysis](https://nbviewer.jupyter.org/github/pinellolab/singlecellvr/blob/master/examples/seurat_integration.ipynb?flush_cache=true). *(Make sure set `n.components = 3` in `RunUMAP(immune.combined, reduction = "pca", dims = 1:30, n.components = 3)`)*

### Velocity:
To get single cell velocity report for scvelo:
``` bash
scvr -t velocity -f pancrease_velocity_updated.h5ad -a clusters
```
* To generate the `pancrease_velocity_updated.h5ad`, check out examples/velocity_3d.ipynb and examples/velocity_3d_scale.ipynb.

### STREAM:  
To get single cell VR report for STREAM : 
```bash
scvr -f ./stream_result/stream_nestorowa16.pkl -t stream -a annotations.txt -g genes.txt -o stream_report
```
* Input files can be found [here](https://www.dropbox.com/sh/fg84hfdeihielun/AACRcmuAIg9RMU30ChgAZevza?dl=0) 
* To generate the `stream_nestorowa16.pkl`, check out [STREAM analysis](https://nbviewer.jupyter.org/github/pinellolab/singlecellvr/blob/master/examples/stream_nestorowa16.ipynb?flush_cache=true).

Or use STREAM package, e.g.:
```python
import stream as st
st.save_vr_report(adata,
                  ann_list=['label','kmeans','branch_id_alias','S4_pseudotime'],
                  gene_list=['Gata1','Car2','Epx','Mfsd2b','Mpo','Emb','Flt3','Dntt'],
                  file_name='stream_report')
```
