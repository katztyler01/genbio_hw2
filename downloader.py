import tiledb
import numpy as np
import pandas as pd
import scanpy as sc
import cellxgene_census
import os
from typing import List


class Downloader:
    def __init__(
        self,
        output_dir: str,
        tissue: str,
        census_version: str = "2023-07-25",
        diseases: List[str] = ["Alzheimer disease", "dementia"],
        cell_types: List[str] = [
            "microglial cell",
            "oligodendrocyte",
            "oligodendrocyte precursor cell",
        ],
    ) -> None:
        self.output_dir = output_dir
        self.tissue = tissue
        self.diseases = diseases
        self.cell_types = cell_types
        self.census_version = census_version

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def filter_obs(self, obs) -> pd.DataFrame:
        print(
            f"Filtering for cell types: {self.cell_types} and diseases: {self.diseases} and subsampling normal cells for each cell type..."
        )
        disease_cells = obs[
            obs["disease"].isin(self.diseases) & obs["cell_type"].isin(self.cell_types)
        ]

        normal_cells = obs[
            (obs["disease"] == "normal") & obs["cell_type"].isin(self.cell_types)
        ]

        counts = pd.crosstab(disease_cells["cell_type"], disease_cells["disease"])

        max_counts = {}
        for cell_type in self.cell_types:
            if cell_type in counts.index:
                max_counts[cell_type] = counts.loc[cell_type].max()
            else:
                max_counts[cell_type] = 0

        sampled_cells = []
        sampled_cells.append(disease_cells)

        for cell_type in self.cell_types:
            normal_ct_cells = normal_cells[normal_cells["cell_type"] == cell_type]
            target_count = max_counts[cell_type]

            if len(normal_ct_cells) <= target_count:
                sampled_normal = normal_ct_cells
                print(
                    f"WARNING: Only {len(normal_ct_cells)} normal {cell_type} cells available, taking all"
                )
            else:
                sampled_idx = np.random.choice(
                    normal_ct_cells.index, size=target_count, replace=False
                )
                sampled_normal = normal_ct_cells.loc[sampled_idx]

            sampled_cells.append(sampled_normal)

        final_dataset = pd.concat(sampled_cells)
        return final_dataset

    def get_anndata(self):
        print(f"Opening census version {self.census_version}...")
        census = cellxgene_census.open_soma(census_version="2023-07-25")

        print(f"Getting {self.tissue} metadata...")
        obs = cellxgene_census.get_obs(
            census,
            "homo_sapiens",
            value_filter=f"tissue_general == '{self.tissue}' and is_primary_data == True",
        )

        filtered_obs = self.filter_obs(obs)
        print("Getting metadata on datasets...")
        census_datasets = (
            census["census_info"]["datasets"]
            .read(
                column_names=[
                    "collection_name",
                    "dataset_title",
                    "dataset_id",
                    "soma_joinid",
                ]
            )
            .concat()
            .to_pandas()
        )
        census_datasets = census_datasets.set_index("dataset_id")

        print(f"Getting metadata for datasets contributing to {self.tissue}...")
        dataset_cell_counts = pd.DataFrame(filtered_obs[["dataset_id"]].value_counts())
        dataset_cell_counts = dataset_cell_counts.rename(columns={0: "cell_counts"})
        dataset_cell_counts = dataset_cell_counts.merge(
            census_datasets, on="dataset_id"
        )

        print("Loading gene metadata...")
        var = cellxgene_census.get_var(census, "homo_sapiens")

        print(f"Getting genes measured in all {self.tissue} datasets")
        presence_matrix = cellxgene_census.get_presence_matrix(
            census, "Homo sapiens", "RNA"
        )
        presence_matrix = presence_matrix[dataset_cell_counts.soma_joinid, :]
        genes_measured = presence_matrix.sum(axis=1).A1
        dataset_cell_counts["genes_measured"] = genes_measured

        var_somaid = np.nonzero(
            presence_matrix.sum(axis=0).A1 == presence_matrix.shape[0]
        )[0].tolist()
        var = var.query(f"soma_joinid in {var_somaid}")
        print(f"Found {var.shape[0]} genes measured in filtered {self.tissue} datasets")

        print(
            f"Downloading {filtered_obs.shape[0]} cells from {obs.shape[0]} total unique cells..."
        )
        subsampled_ids = filtered_obs["soma_joinid"].tolist()

        gene_ids = var["soma_joinid"].to_numpy()
        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_coords=subsampled_ids,
            var_coords=gene_ids,
        )

        adata.var_names = adata.var["feature_name"]

        output_path = os.path.join(
            self.output_dir, f"{self.tissue}_adata_subsampled.h5ad"
        )

        print(f"Saving annData to {output_path}")
        adata.write(output_path)
        census.close()
        return
