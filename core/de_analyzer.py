import scanpy as sc
import diffxpy.api as de
import numpy as np
import pandas as pd
import os
import scanpy as sc
import scvi
import scipy
from scipy.sparse import csr_matrix
import torch
from typing import Optional, List, Dict, Any
import dask


class DEAnalyzer:
    def __init__(
        self,
        output_dir: str,
        adata_path: str,
    ) -> None:
        self.output_dir = output_dir
        self.adata_path = adata_path

        self.adata = sc.read_h5ad(self.adata_path)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def train_scvi(self):
        scvi.data.poisson_gene_selection(self.adata)
        self.adata = self.adata[
            :, self.adata.var["highly_variable"]
        ]  # focus on selected genes
        self.adata.layers["counts"] = self.adata.X.copy().tocsr()
        scvi.model.SCVI.setup_anndata(
            self.adata, layer="counts", batch_key="dataset_id"
        )
        model = scvi.model.SCVI(self.adata, gene_likelihood="nb")
        model.train(
            check_val_every_n_epoch=1,
            max_epochs=100,
            early_stopping=True,
            early_stopping_patience=20,
            early_stopping_monitor="elbo_validation",
        )
        rand_int = np.random.randint(0, 100000)
        model.save(f"data/model_{str(rand_int)}", save_anndata=True)
        return model

    def scvi_de(self, save_dir: Optional[str] = None):
        if save_dir:
            model = scvi.model.SCVI.load(save_dir)
        else:
            model = self.train_scvi()

        for cell_type in self.adata.obs["cell_type"].unique():
            for condition in self.adata.obs["disease"].unique():
                if condition == "normal":
                    continue
                cell_idx1 = (self.adata.obs["cell_type"] == cell_type) & (
                    self.adata.obs["disease"] == condition
                )
                cell_idx2 = (self.adata.obs["cell_type"] == cell_type) & (
                    self.adata.obs["disease"] == "normal"
                )
                de_change = model.differential_expression(
                    idx1=cell_idx1, idx2=cell_idx2
                )
                path = os.path.join(
                    self.output_dir, f"DE_{cell_type}_{condition}.parquet"
                )
                de_change.to_parquet(path)

    def run_de(self):
        control_condition = "normal"
        cell_types = self.adata.obs["cell_type"].unique()
        diseases = [
            d for d in self.adata.obs["disease"].unique() if d != control_condition
        ]
        de_results = {}

        for cell_type in cell_types:
            cell_adata = self.adata[self.adata.obs["cell_type"] == cell_type].copy()
            control_mask = cell_adata.obs["disease"] == control_condition
            if sum(control_mask) < 3:
                print(
                    f"Skipping {cell_type} - not enough control samples ({sum(control_mask)})"
                )
                continue

            for disease in diseases:
                print(f"Running DE for {cell_type}, {disease} vs normal")
                disease_mask = cell_adata.obs["disease"] == disease
                if sum(disease_mask) < 3:
                    print(
                        f"Skipping {cell_type}/{disease} - not enough disease samples ({sum(disease_mask)})"
                    )
                    continue

                comp_adata = cell_adata[control_mask | disease_mask].copy()
                comp_adata.obs["condition"] = (
                    comp_adata.obs["disease"] == disease
                ).astype(int)
                group_label = f"{cell_type}_{disease}_vs_normal"
                print(
                    f"Processing {group_label}: {sum(disease_mask)} disease cells vs {sum(control_mask)} normal cells"
                )

                with dask.config.set(**{"array.slicing.split_large_chunks": True}):
                    test = de.test.wald(
                        data=comp_adata,
                        formula_loc="~ 1 + condition",
                        factor_loc_totest="condition",
                    )

                result = test.summary()
                result["cell_type"] = cell_type
                result["disease"] = disease
                result["comparison"] = f"{disease}_vs_{control_condition}"
                de_results[group_label] = result

                output_path = os.path.join(self.output_dir, f"DE_{group_label}.csv")
                result.to_csv(output_path, index=False)

        return de_results
