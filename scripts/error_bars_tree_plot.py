import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import networkx as nx
from plot_utils import setup_plot, setup_axes
import json
import sys


def main(df_file, outfile, nx_json):
    fig = plot_sample_error_bars_tree_figure(df_file, nx_json)
    fig.savefig(outfile, dpi=300, bbox_inches="tight")


def plot_sample_error_bars_tree_figure(df_file, nx_json, figsize=None):
    if figsize is None:
        figsize = (10, 10)
    df = pd.read_table(df_file)
    setup_plot()
    fig = plt.figure(figsize=figsize, layout="constrained")
    gs = gridspec.GridSpec(1, 2, figure=fig)
    error_bars_ax = fig.add_subplot(gs[0, 0])
    tree_ax = fig.add_subplot(gs[0, 1])
    plot_errorbars(df, error_bars_ax)
    plot_clone_tree(nx_json, tree_ax)
    fig.suptitle("Sample: {}".format(list(df["sample_id"].unique())[0]))
    return fig


def plot_clone_tree(nx_json, tree_ax):
    tree_ax.axis('off')
    with open(nx_json) as f:
        nx_json = json.load(f)
    nx_tree = nx.node_link_graph(nx_json, edges="edges")
    pos = nx.nx_agraph.pygraphviz_layout(nx_tree, prog="dot")
    node_sizes = [(nx_tree.nodes[n]["width"] * 72) ** 2 for n in nx_tree.nodes]
    node_colours = nx.get_node_attributes(nx_tree, "fillcolor")
    node_border_colours = nx.get_node_attributes(nx_tree, "color")
    nx.draw_networkx(
        nx_tree,
        pos=pos,
        ax=tree_ax,
        with_labels=False,
        node_size=node_sizes,
        node_color=node_colours.values(),
        linewidths=0.75,
        edgecolors=node_border_colours.values(),
    )


def plot_errorbars(df, error_bars_ax):
    df = df.loc[~df["clone_id"].str.startswith("ancestral_")]
    df = df.sort_values("clone_tree_order_idx", ascending=False)
    error_bars_ax = setup_axes(error_bars_ax)
    for i in df.index:
        x_val = df.loc[i, "mean"]
        xerr_minus = x_val - df.loc[i, "hdi_lower"]
        xerr_plus = df.loc[i, "hdi_upper"] - x_val
        xerr = [[xerr_minus], [xerr_plus]]
        error_bars_ax.errorbar(
            y=df.loc[i, "clone_id"],
            x=df.loc[i, "mean"],
            xerr=xerr,
            fmt='o',
            c=df.loc[i, "clone_colour"],
        )
    error_bars_ax.set_xlabel("Clonal Prevalence")
    error_bars_ax.set_ylabel("Clone ID")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(
        snakemake.input["rho_summary_table"],
        snakemake.output["error_bars_tree_plot"],
        snakemake.input["prev_sized_nodes_tree_json"],
    )
