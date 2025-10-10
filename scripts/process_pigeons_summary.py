import pandas as pd
import sys


def main(in_file, out_file, model_name, sample_id):
    df = pd.read_csv(in_file)

    df = df[["round", "stepping_stone"]]

    df = df.rename(columns={"stepping_stone": "evidence_estimation"})

    round_to_keep = df["round"].max()

    df = df.loc[df["round"] == round_to_keep]

    df["model_name"] = model_name

    df["sample_id"] = sample_id

    df.to_csv(out_file, index=False, sep="\t")


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(
        snakemake.input["raw_pigeons_summary"],
        snakemake.output["pigeons_summary"],
        snakemake.params["model_used"],
        snakemake.params["sample_id"],
    )
