import pandas as pd
import xarray as xr
import arviz as az
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
from plot_utils import setup_plot, setup_axes
import sys
from matplotlib.patches import Patch


# sort_chroms adapted from: https://github.com/Roth-Lab/hapclone-smk/blob/main/scripts/plot_clone_pseudobulk.py
def sort_chroms(chroms):
    numeric = []
    string = []

    if chroms[0].startswith("chr"):
        chr_prefix = True
    else:
        chr_prefix = False

    for c in chroms:
        if chr_prefix:
            c = c.replace("chr", "")
        try:
            numeric.append(int(c))
        except ValueError:
            string.append(c)

    chroms = [str(x) for x in sorted(numeric)] + list(sorted(string))

    if chr_prefix:
        chroms = ["chr{}".format(x) for x in chroms]

    return chroms


def main(mu_summary_file, p_summary_file, ctDNA_file, outfile, sample_id):

    fig = plot_sample_fits(ctDNA_file, mu_summary_file, p_summary_file, sample_id)
    fig.savefig(outfile, dpi=300, bbox_inches="tight")


def plot_sample_fits(ctDNA_file, mu_summary_file, p_summary_file, sample_id):
    setup_plot()
    fig = plt.figure(figsize=(16, 8.5))
    grid = fig.add_gridspec(2, 1, hspace=0.1)
    ctdna_df = pd.read_table(ctDNA_file)
    # ctdna_df["baf"] = ctdna_df["b"] / ctdna_df["d"]
    ctdna_df["baf"] = ctdna_df["a"] / ctdna_df["d"]
    chroms = sort_chroms(ctdna_df["chrom"].unique())
    chroms_size = ctdna_df["chrom"].value_counts()
    width_ratios = [chroms_size[x] for x in chroms]
    merged_mu = build_merged_param_summary_df(ctdna_df, mu_summary_file)
    merged_p = build_merged_param_summary_df(ctdna_df, p_summary_file)

    plot_vals = [("rdr_cor", "mean", "RDR", merged_mu, True), ("baf", "mean", "BAF", merged_p, False)]
    for i, v in enumerate(plot_vals):
        sub_grid = grid[i].subgridspec(1, len(chroms), width_ratios=width_ratios, wspace=0.05)

        plot_fit(
            v[3],
            chroms,
            fig,
            sub_grid,
            title=v[2],
            y_col=v[0],
            y_col_fit=v[1],
            plot_legend=v[4],
        )
    fig.suptitle("Sample: {}".format(sample_id))
    fig.supxlabel("Chromosome")
    grid.tight_layout(fig)
    return fig


def build_merged_param_summary_df(ctdna_df, param_summary_file):
    param_summary_df = pd.read_table(param_summary_file)
    param_summary_df = pd.merge(param_summary_df, ctdna_df, on="bin")
    return param_summary_df


def plot_fit(df, chroms, fig, grid, title=None, y_col="rdr_cor", y_col_fit="mean", plot_legend=False):
    y_max = max(df[y_col].max(), df["hdi_upper"].max())
    y_min = min(df[y_col].min(), df["hdi_lower"].min())

    weeyums_blue = "#1868DB"

    papaya_orng = "#F47600"

    grouped = df.groupby("chrom")

    last_chrom = chroms[-1]

    for i, chrom in enumerate(chroms):
        chrom_df = grouped.get_group(chrom)

        chrom_df = chrom_df.sort_values(by=["start"])
        num_bins = chrom_df.shape[0]
        chrom_df["idx"] = np.arange(num_bins)
        ax = fig.add_subplot(grid[0, i])
        setup_axes(ax)

        ax.scatter(
            chrom_df["idx"],
            chrom_df[y_col],
            c=papaya_orng,
            s=1,
            alpha=0.5,
        )
        ax.scatter(
            chrom_df["idx"],
            chrom_df[y_col_fit],
            c=weeyums_blue,
            s=1,
            alpha=0.5,
        )
        ax.fill_between(
            chrom_df["idx"],
            y1=chrom_df["hdi_lower"],
            y2=chrom_df["hdi_upper"],
            color=weeyums_blue,
            alpha=0.2,
        )
        sb.despine(ax=ax, offset=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.xaxis.grid(False)
        if i != 0:
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis='y', labelleft=False, left=False)
            # ax.set_yticks([])
            # ax.set_yticklabels([])
        else:
            ax.tick_params(axis="x", which="major", labelsize=12)

        ax.set_xticks([num_bins / 2])
        ax.set_xticklabels([chrom.replace("chr", "")], fontsize=12)

        ax.set_ylim(y_min, y_max)

        if plot_legend and chrom == last_chrom:
            handles = [Patch(facecolor=papaya_orng), Patch(facecolor=weeyums_blue)]
            labels = ("Data", "Model")
            ax.legend(handles, labels, loc="upper right", bbox_to_anchor=(1, 1.15))

    if title is not None:
        ax = fig.add_subplot(grid[:])
        ax.axis("off")
        ax.set_title(title)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr
    main(
        snakemake.input["mu_table"],
        snakemake.input["p_table"],
        snakemake.input["ctDNA_file"],
        snakemake.output["fit_plot"],
        snakemake.params["sample_id"],
    )
