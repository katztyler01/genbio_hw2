import tiledb
import numpy as np
import pandas as pd
import scanpy as sc
import cellxgene_census
import os


class Downloader:
    def __init__(self,
            output_dir: str,
            tissue: str,
            census_version: str = "2023-07-25",
            n_cells: int = 100000
                 ) -> None:
        self.output_dir = output_dir
        self.tissue = tissue
        self.n_cells = n_cells
        self.census_version = census_version
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    def get_anndata(self):
        print(f"Opening census version {self.census_version}...")
        census = cellxgene_census.open_soma(census_version="2023-07-25")
        
        print(f"Getting {self.tissue} metadata...")
        obs = cellxgene_census.get_obs(
            census, "homo_sapiens", value_filter=f"tissue_general == '{self.tissue}' and is_primary_data == True"
        )
        
        print("Getting metadata on datasets...")
        census_datasets = (
            census["census_info"]["datasets"]
            .read(column_names=["collection_name", "dataset_title", "dataset_id", "soma_joinid"])
            .concat()
            .to_pandas()
        )
        census_datasets = census_datasets.set_index("dataset_id")
        
        print(f"Getting metadata for datasets contributing to {self.tissue}...")
        dataset_cell_counts = pd.DataFrame(obs[["dataset_id"]].value_counts())
        dataset_cell_counts = dataset_cell_counts.rename(columns={0: "cell_counts"})
        dataset_cell_counts = dataset_cell_counts.merge(census_datasets, on="dataset_id")
        
        print("Loading gene metadata...")
        var = cellxgene_census.get_var(census, "homo_sapiens")
        
        print(f"Getting genes measured in all {self.tissue} datasets")
        presence_matrix = cellxgene_census.get_presence_matrix(census, "Homo sapiens", "RNA")
        presence_matrix = presence_matrix[dataset_cell_counts.soma_joinid, :]
        genes_measured = presence_matrix.sum(axis=1).A1
        dataset_cell_counts["genes_measured"] = genes_measured
        
        var_somaid = np.nonzero(presence_matrix.sum(axis=0).A1 == presence_matrix.shape[0])[0].tolist()
        var = var.query(f"soma_joinid in {var_somaid}")
        print(f"Found {var.shape[0]} genes measured in all {self.tissue} datasets")
        
        print(f"Subsampling {self.n_cells} cells from {obs.shape[0]} total unique cells...")
        subsampled_n = self.n_cells
        subsampled_ids = obs["soma_joinid"].sample(subsampled_n, random_state=1).tolist()
        
        gene_ids = var["soma_joinid"].to_numpy()
        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_coords=subsampled_ids,
            var_coords=gene_ids,
        )

        adata.var_names = adata.var["feature_name"]
        
        output_path = os.path.join(self.output_dir, f'{self.tissue}_adata_{str(self.n_cells)}.h5ad')
        
        print(f"Saving annData to {output_path}")
        adata.write(output_path)
        census.close()
        return 