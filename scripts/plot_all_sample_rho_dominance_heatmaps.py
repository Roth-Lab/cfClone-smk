from rho_dominance_heatmap import plot_rho_dominance_tree_heatmap
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import sys


def main(rho_table_list, tree_newick, sample_ids, out_file):

    with PdfPages(out_file) as pdf:
        for rho_file, sample_id in zip(rho_table_list, sample_ids):
            fig, ax = plt.subplots()
            ax = plot_rho_dominance_tree_heatmap(rho_file, tree_newick, sample_id, ax=ax)
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(
        snakemake.input["rho_tables"],
        snakemake.input["tree_newick"],
        snakemake.params["sample_ids"],
        snakemake.output["rho_dominance_plots"],
    )
