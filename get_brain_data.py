from downloader import Downloader
import argparse


def read_list(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()
    return list(set([line.strip() for line in lines]))


def main(args):
    kwargs = {}
    if args.disease_list is not None:
        diseases = read_list(args.disease_list)
        kwargs["diseases"] = diseases

    if args.cell_type_list is not None:
        cell_types = read_list(args.cell_type_list)
        kwargs["cell_types"] = cell_types

    if args.census_version is not None:
        kwargs["census_version"] = args.census_version

    downloader = Downloader(output_dir=args.output_dir, tissue=args.tissue, **kwargs)
    downloader.get_anndata()


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Script to download data from cellxgene-census")
    parser.add_argument("--output_dir", help="Path to dave directory", default="data")
    parser.add_argument("--tissue", help="Tissue to download", default="brain")
    parser.add_argument(
        "--census_version",
        required=False,
        help="Census version to download",
        default=None,
    )
    parser.add_argument(
        "--disease_list",
        required=False,
        help="Path to list of diseases to subsample",
        default=None,
    )
    parser.add_argument(
        "--cell_type_list",
        required=False,
        help="Path to list of cell types to subsample",
        default=None,
    )
    args = parser.parse_args()
    main(args)
