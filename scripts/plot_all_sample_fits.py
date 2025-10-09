from plot_cfClone_fits import plot_sample_fits
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import sys


def main(sample_mu_list, sample_p_list, sample_ctDNA_list, sample_id_list, out_file):

    with PdfPages(out_file) as pdf:
        for (
            mu_file,
            p_file,
            ctDNA_file,
            sample_id,
        ) in zip(sample_mu_list, sample_p_list, sample_ctDNA_list, sample_id_list):
            fig = plot_sample_fits(ctDNA_file, mu_file, p_file, sample_id)
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(
        snakemake.input["mu_tables"],
        snakemake.input["p_tables"],
        snakemake.input["ctDNA_files"],
        snakemake.params["sample_ids"],
        snakemake.output["fit_plot"],
    )
