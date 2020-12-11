import utils
import os
import scanpy as sc
import scprep
import tempfile

@utils.loader
def load_from_figshare(URL, filename, test=False):
    
    """
    
    Parameters:
    -----------
    URL
        use a URL from figshare.
        Example from the Nestorowa 2016 dataset: "https://ndownloader.figshare.com/files/25555751"
    
    filename
        the name of the .h5ad file to be saved in a temporary directory. No need to add the ".h5ad" suffix.
        Example: "nestorowa"
    test
        if True, will downsample 1/500th of the data.
        default: False
    
    """
    
    if filename[:-4] != ".h5ad"
        filename = filename + ".h5ad"    
    
    if test:
        adata = load_data(test=False)
        adata = adata[:,:500].copy()
        utils.filter_genes_cells(adata)
        sc.pp.subsample(adata, n_obs=500)
        utils.filter_genes_cells(adata)
        
        return adata
    
    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, filename)
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)
            utils.filter_genes_cells(adata)
            
        return adata
