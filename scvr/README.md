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
Usage: scvr [-h] -f FILE -t {paga,seurat,stream} [-a ANNOTATIONS] [-g GENES] [-o OUTPUT]

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
* To generate the `scanpy_10xpbmc.h5ad`, check out [Scanpy analysis](https://nbviewer.jupyter.org/github/pinellolab/singlecellvr/blob/master/examples/paga3d_paul15.ipynb?flush_cache=true). *(Make sure set `n_components=3` in `sc.tl.umap(adata,n_components=3)`)*


### PAGA:  

To get single cell VR report for PAGA :  
```bash
scvr -f ./paga_result/paga3d_paul15.h5ad -t paga -a annotations.txt -g genes.txt -o paga_report
```

* Input files can be found [here](https://www.dropbox.com/sh/03zpxs9zv7yusi1/AADKVSU8Il1JcjA7lfHjmRpSa?dl=0) 
* To generate the `paga3d_paul15.h5ad`, check out [PAGA analysis](https://nbviewer.jupyter.org/github/pinellolab/singlecellvr/blob/master/examples/paga3d_paul15.ipynb?flush_cache=true). *(Make sure set `n_components=3` in `sc.tl.umap(adata,n_components=3)`)*

### Seurat:  
To get single cell VR report for Seurat :  
```bash
scvr -f ./seurat_result/seurat3d_10xpbmc.loom -t seurat -a annotations.txt -g genes.txt -o seurat_report
```
* Input files can be found [here](https://www.dropbox.com/sh/tpk4qfm5qsjpffn/AADmKmyDx7rhzKBOpIlAgMEUa?dl=0) 
* To generate the `seurat3d_10xpbmc.loom`, check out [Seurat analysis](https://nbviewer.jupyter.org/github/pinellolab/singlecellvr/blob/master/examples/seurat3d_10xpbmc.ipynb?flush_cache=true). *(Make sure set `n.components = 3` in `pbmc <- RunUMAP(pbmc, dims = 1:10, n.components = 3)`)*

### STREAM:  
To get single cell VR report for STREAM : 
```bash
scvr -f ./stream_result/stream_nestorowa16.pkl -t stream -a annotations.txt -g genes.txt -o stream_report
```
* Input files can be found [here](https://www.dropbox.com/sh/fg84hfdeihielun/AACRcmuAIg9RMU30ChgAZevza?dl=0) 
* To generate the `stream_nestorowa16.pkl`, check out [STREAM analysis](https://nbviewer.jupyter.org/github/pinellolab/singlecellvr/blob/master/examples/stream_nestorowa16.ipynb?flush_cache=true).