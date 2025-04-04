from downloader import Downloader
import argparse

def main(args):
    downloader = Downloader(
        output_dir=args.output_dir, 
        tissue=args.tissue,
        census_version=args.census_version,
        n_cells=args.n_cells
    )
    downloader.get_anndata()


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Script to download data from cellxgene-census")
    parser.add_argument("--output_dir", help="Path to dave directory", default="data")
    parser.add_argument("--tissue", help="Tissue to download", default="brain")
    parser.add_argument("--census_version", required=False, help="Census version to download", default="2023-07-25")
    parser.add_argument("--n_cells", required=False, help="Number of cells to subsample", default=100000)
    args = parser.parse_args()
    main(args)