import pandas as pd
from neo4j import GraphDatabase
import numpy as np
import dotenv
import os
import argparse

# NOTE: code adapted from https://neo4j.com/docs/python-manual/current/


def load_data(parquet_file):
    df = pd.read_parquet(parquet_file)
    gene_nodes = []
    disease_nodes = []
    cell_type_edges = []

    for _, row in df.iterrows():
        row_dict = {}
        for k, v in row.items():
            if isinstance(v, (list, np.ndarray)):
                row_dict[k] = v
            else:
                row_dict[k] = None if pd.isna(v) else v

        gene = {
            "symbol": row_dict["feature_name"],
            "entrez": row_dict["entrez"],
            "ensembl": row_dict["ensembl"],
        }
        gene_nodes.append(gene)
        disease = {
            "name": row_dict["disease"],
            "mondo": row_dict["mondo"],
            "doid": row_dict["doid"],
            "mesh": row_dict["mesh"],
        }
        disease_nodes.append(disease)
        cell_type_edge = {
            "gene_symbol": row_dict["feature_name"],
            "disease_name": row_dict["disease"],
            "cell_type": row_dict["cell_type"],
            "lfc_mean": row_dict["lfc_mean"],
        }
        cell_type_edges.append(cell_type_edge)
    return gene_nodes, disease_nodes, cell_type_edges


def main(args):
    load_status = dotenv.load_dotenv(args.env_file)
    if load_status is False:
        raise RuntimeError("Environment variables not loaded.")

    URI = os.getenv("NEO4J_URI")
    AUTH = (os.getenv("NEO4J_USERNAME"), os.getenv("NEO4J_PASSWORD"))

    with GraphDatabase.driver(URI, auth=AUTH) as driver:
        driver.verify_connectivity()
        print("Connection established.")

        gene_nodes, disease_nodes, cell_type_edges = load_data(args.de_parquet)
        print("Data loaded from parquet file.")

        for gene in gene_nodes:
            driver.execute_query(
                """
                MERGE (g:Gene {symbol: $gene.symbol})
                ON CREATE SET 
                    g.entrez = $gene.entrez,
                    g.ensembl = $gene.ensembl
                """,
                gene=gene,
                database_="neo4j",
            )
        for disease in disease_nodes:
            driver.execute_query(
                """
                MERGE (d:Disease {name: $disease.name})
                ON CREATE SET 
                    d.mondo = $disease.mondo,
                    d.doid = $disease.doid,
                    d.mesh = $disease.mesh
                """,
                disease=disease,
                database_="neo4j",
            )
        for edge in cell_type_edges:
            cell_type = edge["cell_type"]
            clean_cell_type = "".join(
                c if c.isalnum() else "_" for c in cell_type
            ).strip("_")

            driver.execute_query(
                """
                MATCH (g:Gene {symbol: $rel.gene_symbol})
                MATCH (d:Disease {name: $rel.disease_name})
                CALL apoc.create.relationship(g, $rel_type, {lfc_mean: $rel.lfc_mean}, d)
                YIELD rel
                RETURN rel
                """,
                rel=edge,
                rel_type=f"DE_IN_{clean_cell_type}",
                database_="neo4j",
            )
        print("Data loaded into Neo4j database.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Load data into Neo4j database.")
    parser.add_argument(
        "--de_parquet", type=str, help="Path to the input parquet file."
    )
    parser.add_argument(
        "--env_file", type=str, help="Path to the environment variables file."
    )
    args = parser.parse_args()

    main(args)

# NOTE: To visualize top 5 genes for each disease in Neo4j browser, run this query:
"""
MATCH (g:Gene)-[r]->(d:Disease)
WHERE type(r) STARTS WITH 'DE_IN_'
WITH g, d, r, type(r) as relationship_type, abs(r.lfc_mean) as abs_lfc
WITH g, d, r, relationship_type, abs_lfc,
     substring(relationship_type, 6) as cell_type  // Remove 'DE_IN_' prefix

WITH g, d, collect({rel: r, cell_type: cell_type, abs_lfc: abs_lfc}) as edges,
     max(abs_lfc) as max_abs_lfc

ORDER BY d.name, max_abs_lfc DESC
WITH d, collect({gene: g, edges: edges, max_abs_lfc: max_abs_lfc})[0..5] as top_genes
UNWIND top_genes as gene_data
UNWIND gene_data.edges as edge_data

RETURN gene_data.gene as gene, 
       edge_data.rel as relationship, 
       d as disease, 
       edge_data.cell_type as cell_type,
       edge_data.abs_lfc as abs_lfc
"""
