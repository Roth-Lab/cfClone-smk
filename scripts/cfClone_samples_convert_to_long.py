import pandas as pd
import sys


def build_rho_cluster_df_from_raw_samples(infile):
    df = pd.read_csv(infile)

    df = pd.wide_to_long(df, stubnames="rho", sep="_", i=["iteration", "chain"], j="clone_id", suffix='(\\d+|normal)')

    df.reset_index(inplace=True)

    keep_cols = [col for col in df.columns if not col.startswith("mu.") and not col.startswith("p.")]

    df = df[keep_cols]
    return df


def main(infile, out_clonal_rho):
    clonal_df = build_rho_cluster_df_from_raw_samples(infile)

    clonal_df.to_csv(out_clonal_rho, sep="\t", index=False)


if __name__ == '__main__':
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(snakemake.input["cfClone_samples"], snakemake.output["cfClone_clonal_rho_table"])
