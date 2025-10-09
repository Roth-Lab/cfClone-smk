import pandas as pd
import xarray as xr
import arviz as az
from create_rho_summary import define_hdi_upper_and_lower_cols
import sys


def build_mu_p_cluster_df_from_raw_samples(infile):
    df = pd.read_csv(infile)

    df = pd.wide_to_long(df, stubnames=["mu", "p"], sep=".", i=["iteration", "chain"], j="bin")

    df.reset_index(inplace=True)
    return df


def build_bin_col(df_summary):
    df_summary = df_summary.reset_index(names=["parameters"])
    df_summary[["parameter_name", "bin"]] = df_summary["parameters"].str.split('[', expand=True)
    df_summary["bin"] = df_summary["bin"].str.removesuffix(']')
    df_summary["bin"] = pd.to_numeric(df_summary["bin"])
    return df_summary


def build_arviz_summary_df(df, varname):
    df = df.rename(columns={"iteration": "draw"})
    df = df.set_index(["chain", "draw", "bin"])
    xdata = xr.Dataset.from_dataframe(df)
    az_dataset = az.InferenceData(posterior=xdata)
    df_summary = az.summary(az_dataset, var_names=[varname], kind="stats", hdi_prob=0.95, round_to="none")
    return df_summary


def build_parameter_summary_df(df, varname):
    param_summary_df = build_arviz_summary_df(df, varname)
    param_summary_df = build_bin_col(param_summary_df)
    define_hdi_upper_and_lower_cols(param_summary_df)
    return param_summary_df


def main(infile, mu_summary_out, p_summary_out):
    df = build_mu_p_cluster_df_from_raw_samples(infile)

    mu_df = build_parameter_summary_df(df, "mu")
    p_df = build_parameter_summary_df(df, "p")

    mu_df.to_csv(mu_summary_out, sep="\t", index=False)
    p_df.to_csv(p_summary_out, sep="\t", index=False)


if __name__ == '__main__':
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(snakemake.input["cfClone_samples"], snakemake.output["mu_table"], snakemake.output["p_table"])
