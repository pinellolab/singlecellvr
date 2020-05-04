SingleCellVR Preprocess:  
Prepare your data for the visualization on single cell VR website <https://singlecellvr.com/>

Installation
------------
Install and update using pip:  
`pip install scvr-prep`

Usage
-----
`$ scvr_prep --help`

```
- Single cell VR preprocessing -
Version 1.1

usage: scvr_prep [-h] -f FILE [-g GENES] [-a ANNOTATIONS] -t
                 {paga,seurat,stream} [-o OUTPUT]

scvr_prep Parameters

required arguments:
  -f FILE, --filename FILE
                        Analysis result file name (default: None)
  -t {paga,seurat,stream}, --toolname {paga,seurat,stream}
                        Tool name used to generate the analysis result.
                        (default: None)
                        
optional arguments:
  
  -g GENES, --genes GENES
                        Gene list file name. It contains the genes to
                        visualize in one column (default: None)
  -a ANNOTATIONS, --annotations ANNOTATIONS
                        Annotation file name. It contains cell annotations
                        used to color cells (default: None)
  -o OUTPUT, --output OUTPUT
                        Output folder (default: vr_report) (default:
                        vr_report)
  -h, --help            show this help message and exit                        
```


Examples:
---------
* PAGA:  
`scvr_prep -f ./paga_result/paga3d_paul15.h5ad -t paga -a annotations.txt -g genes.txt -o paga_report`
> To generate the `paga3d_paul15.h5ad`, check out PAGA analysis script: 

* Seurat:  
`scvr_prep -f ./seurat_result/seurat3d_10xpbmc.loom -a annotations.txt -g genes.txt -o seurat_report -t seurat`

* STREAM:  
`scvr_prep -f ./stream_result/stream_result.pkl -g genes.txt -o stream_report -t stream`