import pandas as pd
import pathlib


def main(args):
    df = pd.read_csv(args.in_file, sep="\t")

    out_dir = pathlib.Path(args.out_dir)

    out_dir.mkdir(exist_ok=True, parents=True)

    for clone in df["clone"].unique():
        open(out_dir.joinpath(f"clone_{clone}.txt"), "w").close()

    for r in ["full", "normal"]:
        open(out_dir.joinpath(f"{r}.txt"), "w").close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-dir", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
