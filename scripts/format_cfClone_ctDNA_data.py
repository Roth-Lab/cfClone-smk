import pandas as pd
import sys


def main(in_file, out_file, in_cn_profile, out_cn_profile):
    ctdna_df = get_ctdna_df(in_file)

    add_bin_name_col(ctdna_df)

    in_cn_df = pd.read_table(in_cn_profile)

    add_bin_name_col(in_cn_df)

    bin_names_intersect = set(ctdna_df["bin_name"].unique()).intersection(in_cn_df["bin_name"].unique())

    in_cn_df = in_cn_df.loc[in_cn_df["bin_name"].isin(bin_names_intersect)]
    ctdna_df = ctdna_df.loc[ctdna_df["bin_name"].isin(bin_names_intersect)]

    in_cn_df.drop(columns=["bin_name"], inplace=True)
    ctdna_df.drop(columns=["bin_name"], inplace=True)

    ctdna_df["bin"] = range(1, len(ctdna_df) + 1)

    in_cn_df.to_csv(out_cn_profile, index=False, sep="\t")
    ctdna_df.to_csv(out_file, index=False, sep="\t")


def add_bin_name_col(df):
    df["bin_name"] = df["chrom"].astype(str) + "_" + df["start"].astype(str) + "_" + df["end"].astype(str)


def get_ctdna_df(in_file):
    df = pd.read_table(in_file)
    df = df.loc[df["valid"] == True]
    df = df.rename(columns={"allele_0_count": "a", "allele_1_count": "b"})
    df["d"] = df["a"] + df["b"]
    df = df.sort_values(by=["chrom", "start", "end"], ignore_index=True)
    df = df[['chrom', 'start', 'end', 'a', 'b', 'd', 'rdr_cor']].copy()
    return df


if __name__ == '__main__':
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(
        snakemake.input["raw_ctDNA_file"],
        snakemake.output["ctDNA_file"],
        snakemake.input["patient_cn_profiles"],
        snakemake.output["cn_profile"],
    )
