from skbio.tree import TreeNode

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd

from phylo_clustergrid import clustermap_phylo


def main(args):
    df = pd.read_csv(args.in_file, index_col="clone_id", sep="\t")

    df.index = df.index.astype(str)

    tree = TreeNode.read(args.tree_file, convert_underscores=False)

    cmap_greater_than = mcolors.LinearSegmentedColormap.from_list(
        "red_white_blue", [(0, "blue"), (0.5, "white"), (1, "red")]
    )

    fig = clustermap_phylo(
        df,
        cmap=cmap_greater_than,
        cbar_kws={'label': r"$P(\rho_{i} > \rho_{j})$"},
        col_cluster=False,
        row_tree=tree,
        use_branch_lengths=False,
        label=True,
    )

    fig.ax_heatmap.set(
        xlabel='Clone ID: $\\mathbf{i}$', ylabel='Clone ID: $\\mathbf{j}$'
    )

    if args.sample_id is not None:
        fig.ax_col_dendrogram.set(title="Sample: {}".format(args.sample_id))

    plt.tight_layout()

    plt.savefig(args.out_file, dpi=300, bbox_inches='tight')

    plt.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-t", "--tree-file", required=True)

    parser.add_argument("-s", "--sample-id", default=None)

    cli_args = parser.parse_args()

    main(cli_args)
