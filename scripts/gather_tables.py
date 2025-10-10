import pandas as pd
import sys


def main(in_files, out_file):
    df = pd.concat(map(pd.read_table, in_files))
    df.to_csv(out_file, index=False, sep="\t")


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(snakemake.input, snakemake.output["out_table"])
