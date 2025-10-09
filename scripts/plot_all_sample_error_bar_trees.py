from error_bars_tree_plot import plot_sample_error_bars_tree_figure
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import sys


def main(rho_table_list, tree_json_list, out_file):

    with PdfPages(out_file) as pdf:
        for rho_file, tree_json_file in zip(rho_table_list, tree_json_list):
            fig = plot_sample_error_bars_tree_figure(rho_file, tree_json_file)
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(snakemake.input["rho_tables"], snakemake.input["tree_jsons"], snakemake.output["errorbar_tree_plot"])
