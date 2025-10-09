import sys
import pandas as pd


def main(rho_tables, out_file):
    df = pd.concat(map(pd.read_table, rho_tables))
    df.drop(columns=["parameters", "hdi_lower", "hdi_upper", "clone_tree_order_idx", "clone_colour"], inplace=True)
    cols_not_to_rename = {"clone_id", "sample_id", "parameter_name"}
    col_names_update_dict = {
        col: "{}{}".format("clonal_prevalence_", col) for col in df.columns if col not in cols_not_to_rename
    }
    df.rename(columns=col_names_update_dict, inplace=True)
    df.to_csv(out_file, index=False, sep="\t")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(snakemake.input["rho_tables"], snakemake.output["clonal_prevalence_table"])
