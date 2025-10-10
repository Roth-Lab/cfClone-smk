import pandas as pd
import sys


def main(in_file, out_file, clone):
    df = pd.read_table(in_file)
    clones_to_keep = {"normal"}
    clones_to_keep.add(clone)
    df = df.loc[df["clone"].isin(clones_to_keep)]
    df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(snakemake.input["cn_profile"], snakemake.output["cn_profile"], snakemake.params["clone"])
