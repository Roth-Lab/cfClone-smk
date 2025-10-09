import pandas as pd
import sys


def main(bed_file, out_file, add_normal=True):
    df = get_initial_clone_cn_df(bed_file)

    if add_normal:
        df = add_normal_clone(df)

    df.to_csv(out_file, index=False, sep="\t")


def add_normal_clone(df):
    new_df = df[["chrom", "start", "end"]].drop_duplicates(ignore_index=True)
    new_df["clone"] = "normal"
    new_df["cn_a"] = 1
    new_df["cn_b"] = 1
    new_df["cn_t"] = 2
    new_df["cn_baf"] = 1 / 2
    merged_df = pd.concat([df, new_df])
    merged_df.sort_values(by=["chrom", "start", "end", "clone"], inplace=True, ignore_index=True)
    return merged_df


def get_initial_clone_cn_df(bed_file):
    bed_df = pd.read_table(bed_file)
    long_df = pd.wide_to_long(bed_df, stubnames='', i=["chr", "start", "end"], j="clone")
    long_df.reset_index(inplace=True)
    long_df.rename(columns={"chr": "chrom", "": "phased_cn"}, inplace=True)
    spit_cn_cols = ["cn_a", "cn_b"]
    long_df[spit_cn_cols] = long_df["phased_cn"].str.split("|", expand=True)
    long_df.drop(columns=["phased_cn"], inplace=True)
    long_df[spit_cn_cols] = long_df[spit_cn_cols].apply(pd.to_numeric)
    long_df["cn_t"] = long_df["cn_a"] + long_df["cn_b"]
    long_df["cn_baf"] = long_df["cn_b"] / long_df["cn_t"]
    return long_df


if __name__ == '__main__':
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(snakemake.input["bed_file"], snakemake.output["formatted_cn_profiles"])
