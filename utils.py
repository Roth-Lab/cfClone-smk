import pandas as pd
from pathlib import Path


def _get_default_cfclone_settings_dict():
    default_params = {
        "subsampling": 1,
        "n_rounds": 10,
        "num_threads": 10,
        "symmetrization": 0.5,
        "rdr_likelihood": "TDistribution",
        "baf_likelihood": "BetaBinomial",
        "keep_only_normal": False,
    }
    return default_params


class ConfigManager(object):
    def __init__(self, config):
        self.config = config

        self.out_dir = Path(self.config["out_dir"]).resolve()

        self.log_dir = self.out_dir.joinpath("log")

        self.benchmarking_dir = self.out_dir.joinpath("benchmarks")

        self.sample_df = self._load_sample_df()

        self.sample_list = list(self.sample_df["sample_id"].unique())

        self.sample_df.set_index("sample_id", inplace=True)

        self.csv_path = Path("").joinpath("results", "latest", "samples.csv.gz")

        self.raw_pigeons_summary_path = Path("").joinpath("results", "latest", "build", "Pigeons_summary.csv")

        self._load_cfclone_settings()

        self.cfclone_settings = self.config["cfClone"]

        self.clones = self.get_clone_list()

    def get_clone_list(self):
        bed_df = pd.read_table(self.patient_cn_bed_file)
        non_clone_cols = {"chr", "start", "end"}
        clone_list = [col for col in bed_df.columns if col not in non_clone_cols]
        clone_list.append("normal")
        return clone_list

    def _load_cfclone_settings(self):
        default_params = _get_default_cfclone_settings_dict()
        cfclone_settings = self.config.get("cfClone", {})
        default_params.update(cfclone_settings)
        self.config["cfClone"] = default_params

    @staticmethod
    def add_suffix(file_path, suffix_to_add):
        return file_path.with_suffix(file_path.suffix + suffix_to_add)

    @property
    def raw_clone_tree(self):
        return Path(self.config["clone_tree_newick"])

    def get_log_file(self, rule_name, file_name):
        log_path = self.log_dir.joinpath(rule_name, file_name)
        file_name_suffix = log_path.suffix
        if file_name_suffix != ".log":
            log_path = self.add_suffix(log_path, ".log")
        return log_path

    def get_benchmark_file(self, rule_name, file_name):
        log_path = self.benchmarking_dir.joinpath(rule_name, file_name)
        return log_path

    def _load_sample_df(self):
        df = pd.read_table(self.config["sample_tsv"])
        return df

    def get_sample_raw_ctDNA_file(self, wildcards):
        return self.sample_df.loc[wildcards.sample, "path"]

    def get_restricted_model_name(self, wildcards):
        clone = wildcards.clone
        if clone == "normal":
            return "normal_only"
        else:
            return "normal_and_clone_{}".format(clone)

    @property
    def cfclone_threads(self):
        return self.cfclone_settings["num_threads"]

    @property
    def subsampling(self):
        return self.cfclone_settings["subsampling"]

    @property
    def n_rounds(self):
        return self.cfclone_settings["n_rounds"]

    @property
    def symmetrization(self):
        return self.cfclone_settings["symmetrization"]

    @property
    def rdr_likelihood(self):
        return self.cfclone_settings["rdr_likelihood"]

    @property
    def baf_likelihood(self):
        return self.cfclone_settings["baf_likelihood"]

    @property
    def keep_only_normal(self):
        keep_only_normal = self.cfclone_settings["keep_only_normal"]
        if keep_only_normal:
            return "--keep-only-normal"
        else:
            return ""

    @property
    def samples_dir(self):
        return self.out_dir.joinpath("samples")

    @property
    def sample_run_dir(self):
        return self.samples_dir.joinpath("{sample}")

    @property
    def restricted_models_dir(self):
        return self.sample_run_dir.joinpath("restricted_models")

    @property
    def restricted_model_run_dir(self):
        return self.restricted_models_dir.joinpath("{clone}")

    @property
    def restricted_cfclone_sample_results_dir(self):
        return self.restricted_model_run_dir.joinpath("cfclone_runs")

    @property
    def restricted_cfclone_results_file(self):
        return self.restricted_model_run_dir.joinpath("cfclone_results.csv.gz")

    @property
    def raw_restricted_pigeons_summary(self):
        return self.restricted_model_run_dir.joinpath("Pigeons_summary.csv")

    @property
    def restricted_pigeons_summary(self):
        return self.restricted_model_run_dir.joinpath("clone_{clone}_restricted_pigeons_summary.tsv.gz")

    # @property
    # def restricted_sample_ctDNA_file(self):
    #     return self.restricted_model_run_dir.joinpath("inputs", "restricted_cfClone_ctDNA.tsv.gz")

    @property
    def restricted_sample_cn_profile(self):
        return self.restricted_model_run_dir.joinpath("restricted_sample_cn_profiles.tsv.gz")

    @property
    def run_input_dir(self):
        return self.sample_run_dir.joinpath("inputs")

    @property
    def cfclone_sample_results_dir(self):
        return self.sample_run_dir.joinpath("cfclone_runs")

    @property
    def cfclone_results_file(self):
        return self.sample_run_dir.joinpath("cfclone_results.csv.gz")

    @property
    def raw_pigeons_summary(self):
        return self.sample_run_dir.joinpath("Pigeons_summary.csv")

    @property
    def pigeons_summary(self):
        return self.sample_run_dir.joinpath("full_model_pigeons_summary.tsv.gz")

    @property
    def sample_pigeons_summary(self):
        return self.sample_run_dir.joinpath("all_models_pigeons_summary.tsv.gz")

    @property
    def sample_ctDNA_file(self):
        return self.run_input_dir.joinpath("cfClone_ctDNA.tsv.gz")

    @property
    def sample_cn_profile(self):
        return self.run_input_dir.joinpath("sample_cn_profiles.tsv.gz")

    @property
    def sample_processed_output_dir(self):
        return self.sample_run_dir.joinpath("processed_outputs")

    @property
    def cfclone_samples_clone_long(self):
        return self.sample_processed_output_dir.joinpath("cfClone_samples_clone_long.tsv.gz")

    @property
    def cfclone_samples_clonal_rho(self):
        return self.sample_processed_output_dir.joinpath("cfClone_samples_clonal_rho.tsv.gz")

    @property
    def cfClone_p_table(self):
        return self.sample_processed_output_dir.joinpath("cfClone_p_table.tsv.gz")

    @property
    def cfClone_mu_table(self):
        return self.sample_processed_output_dir.joinpath("cfClone_mu_table.tsv.gz")

    @property
    def sample_clonal_rho_summary(self):
        return self.sample_processed_output_dir.joinpath("clonal_rho_summary.tsv.gz")

    @property
    def prev_sized_nodes_tree_json(self):
        return self.sample_processed_output_dir.joinpath("prev_sized_nodes_tree.json")

    @property
    def labelled_nodes_tree_json(self):
        return self.sample_processed_output_dir.joinpath("labelled_nodes_tree.json")

    @property
    def patient_cn_profiles(self):
        return self.out_dir.joinpath("patient_cn_profiles.tsv.gz")

    @property
    def patient_cn_bed_file(self):
        return Path(self.config["clonal_cn_bed"])

    @property
    def clonal_prevalence_table(self):
        return self.out_dir.joinpath("clonal_prevalence_table.tsv.gz")

    @property
    def all_pigeons_summary_table(self):
        return self.out_dir.joinpath("all_samples_pigeons_summary.tsv.gz")

    @property
    def experiment_configuration(self):
        return self.out_dir.joinpath("experiment_configuration.yaml")

    @property
    def plot_dir(self):
        return self.out_dir.joinpath("plots")

    @property
    def fits_plot(self):
        return self.plot_dir.joinpath("fits_plot.pdf")

    @property
    def rho_dominance_plot(self):
        return self.plot_dir.joinpath("rho_dominance_plot.pdf")

    @property
    def gathered_error_bars_tree_plot(self):
        return self.plot_dir.joinpath("errorbars_tree_plot.pdf")

    @property
    def plot_resource_dir(self):
        return self.plot_dir.joinpath("plot_resources")

    @property
    def all_sample_plots_dir(self):
        return self.plot_dir.joinpath("sample_plots")

    @property
    def sample_plots_dir(self):
        return self.all_sample_plots_dir.joinpath("{sample}")

    @property
    def error_bars_tree_plot(self):
        return self.sample_plots_dir.joinpath("error_bars_tree_plot.svg")

    @property
    def rdr_baf_fit_plot(self):
        return self.sample_plots_dir.joinpath("rdr_baf_fit_plot.svg")

    @property
    def sample_rho_dominance_heatmap(self):
        return self.sample_plots_dir.joinpath("rho_dominance_heatmap.png")

    @property
    def prev_sized_nodes_tree_svg(self):
        return self.sample_plots_dir.joinpath("prev_sized_nodes_tree.svg")

    @property
    def labelled_nodes_tree_png(self):
        return self.sample_plots_dir.joinpath("labelled_nodes_tree.png")

    @property
    def labelled_newick(self):
        return self.plot_resource_dir.joinpath("preprocessed_clone_tree.nwk")

    @property
    def colour_dict_json(self):
        return self.plot_resource_dir.joinpath("clonal_colour_dict.json")
