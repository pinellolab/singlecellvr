singlecellVR preprocess: prepare your data for the visualization on single cell VR website https://singlecellvr.com/

`scvr_prep -f ./stream_result/stream_result.pkl -g genes.txt -o stream_report -t stream`

`scvr_prep -f ./paga_result/paga3d_paul15.h5ad -a annotations.txt -g genes.txt -o paga_report -t paga`

`scvr_prep -f ./seurat_result/seurat3d_10xpbmc.loom -a annotations.txt -g genes.txt -o seurat_report -t seurat `