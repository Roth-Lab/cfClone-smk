import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
from skbio.tree import TreeNode
import sys
from plot_utils import setup_plot
import seaborn as sb


def main(clonal_rho_table, tree_heatmap, treefile, sample_id):
    fig_ax = plot_rho_dominance_tree_heatmap(clonal_rho_table, treefile, sample_id)
    plt.tight_layout()
    plt.savefig(tree_heatmap, dpi=300, bbox_inches='tight')
    plt.close()


def plot_rho_dominance_tree_heatmap(clonal_rho_table, treefile, sample_id, ax=None):
    setup_plot()
    df = pd.read_table(clonal_rho_table)
    new_df2 = build_dominance_df(df)
    row_tree = TreeNode.read(treefile, convert_underscores=False)

    leaf_idx_map = {tip.name: i for i, tip in enumerate(row_tree.tips())}

    new_df2["tree_idx"] = leaf_idx_map

    new_df2.sort_values("tree_idx", inplace=True)

    new_df2.drop(columns="tree_idx", inplace=True)

    new_col_order = new_df2.index.tolist()

    new_df2 = new_df2[new_col_order]

    cmap_greater_than = mcolors.LinearSegmentedColormap.from_list(
        "red_white_blue", [(0, "blue"), (0.5, "white"), (1, "red")]
    )

    heat_ax = sb.heatmap(
        new_df2,
        cmap=cmap_greater_than,
        cbar_kws={'label': r"$P(\rho_{i} > \rho_{j})$"},
        ax=ax,
    )

    heat_ax.set(
        xlabel='Clone ID: $\\mathbf{i}$',
        ylabel='Clone ID: $\\mathbf{j}$',
        title="Sample: {}".format(sample_id),
    )

    return heat_ax


def build_dominance_df(df):
    df = df[["iteration", "chain", "clone_id", "rho"]]
    df = df.set_index(["iteration", "chain"])
    clonal_grouped = df.groupby("clone_id")
    new_df = []
    for i_name, i_group in clonal_grouped:
        i_records = {"clone_id": i_name}
        total_len = len(i_group)
        for j_name, j_group in clonal_grouped:
            diff_df = i_group["rho"] > j_group["rho"]
            diff_val = diff_df.sum() / total_len
            i_records[j_name] = diff_val
        new_df.append(i_records)
    new_df = pd.DataFrame(new_df)
    new_df2 = new_df.set_index("clone_id")
    return new_df2


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(
        snakemake.input["clonal_rho_table"],
        snakemake.output["rho_heatmap"],
        snakemake.input["processed_clone_tree"],
        snakemake.params["sample_id"],
    )
