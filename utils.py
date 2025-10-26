from pathlib import Path

import pandas as pd


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

        # self.sample_df = self._load_sample_df()

        # self.sample_list = list(self.sample_df["sample_id"].unique())

        # self.sample_df.set_index("sample_id", inplace=True)

        # self.csv_path = Path("").joinpath("results", "latest", "samples.csv.gz")

        # self.raw_pigeons_summary_path = Path("").joinpath(
        #     "results", "latest", "build", "Pigeons_summary.csv"
        # )

        # self._load_cfclone_settings()

        # self.cfclone_settings = self.config["cfClone"]

        # self.clones = self.get_clone_list()

    # Params
    @property
    def clones(self):
        df = pd.read_csv(self.clone_cn_file, sep="\t")
        return df["clone"].unique()

    @property
    def samples(self):
        return self.config["ctdna_files"].keys()

    # cfClone params
    @property
    def num_chains(self):
        return self.config["num_chains"]

    @property
    def num_chains_vi(self):
        return self.config["num_chains_vi"]

    @property
    def num_rounds(self):
        return self.config["num_rounds"]

    @property
    def num_threads(self):
        return self.config["num_threads"]

    @property
    def post_process_results(self):
        return [
            "dominance_prob",
            "pairwise_ranks",
            "parameter_summaries",
            "prevalence_samples",
            "prevalence_stats",
            "samples",
            "summary",
            "tumour_content",
        ]

    @property
    def run_types(self):
        return ["full", "normal"] + ["clone_{}".format(x) for x in self.clones]

    # Directories
    @property
    def benchmark_dir(self):
        return self.pipeline_dir.joinpath("benchmark")

    @property
    def log_dir(self):
        return self.pipeline_dir.joinpath("log")

    @property
    def out_dir(self):
        return Path(self.config["out_dir"])

    @property
    def pipeline_dir(self):
        return Path(self.config["pipeline_dir"])

    @property
    def sample_out_dir(self):
        return self.out_dir.joinpath("{sample}")

    @property
    def sample_tmp_dir(self):
        return self.tmp_dir.joinpath("{sample}")

    @property
    def tmp_dir(self):
        return self.pipeline_dir.joinpath("tmp")

    # Input files
    @property
    def clone_cn_file(self):
        return Path(self.config["clone_cn_file"]).resolve()

    @property
    def clone_tree_file(self):
        return Path(self.config["clone_tree_newick"]).resolve()

    def get_ctdna_file(self, wc):
        return self._get_ctdna_file(wc.sample)

    # Pipeline files
    @property
    def ancestral_prevalence_template(self):
        return self.sample_out_dir.joinpath("tables", "ancestral_prevalence.tsv.gz")

    @property
    def clone_prevalence_tree_template(self):
        return self.sample_out_dir.joinpath("trees", "prevalence_tree.json")

    @property
    def clone_prevalences_plot(self):
        return self.sample_out_dir.joinpath("plots", "clone_prevalences.pdf")

    @property
    def evidence_template(self):
        return self.sample_out_dir.joinpath("tables", "evidence.tsv")

    @property
    def experiment_configuration(self):
        return self.out_dir.joinpath("config.yaml")

    @property
    def fit_template(self):
        return self.sample_out_dir.joinpath("fit", "{run_type}.h5")

    @property
    def fit_plot_template(self):
        return self.sample_out_dir.joinpath("plots", "fit.pdf")

    @property
    def pairwise_ranks_plot(self):
        return self.sample_out_dir.joinpath("plots", "pairwise_ranks.pdf")

    @property
    def post_process_template(self):
        return self.sample_out_dir.joinpath("tables", "{pp_result}.tsv")

    @property
    def run_type_evidence_template(self):
        return self.sample_tmp_dir.joinpath("evidence", "{run_type}.csv")

    @property
    def pipeline_files(self):
        files = []

        files.append(self.experiment_configuration)

        for s in self.samples:
            files.append(str(self.clone_prevalences_plot).format(sample=s))

            files.append(str(self.evidence_template).format(sample=s))

            files.append(str(self.fit_plot_template).format(sample=s))

            files.append(str(self.pairwise_ranks_plot).format(sample=s))

            for r in self.run_types:
                if r == "full":
                    for p in self.post_process_results:
                        files.append(
                            str(self.post_process_template).format(
                                pp_result=p, sample=s
                            )
                        )

        return files

    # Helper functions
    def get_benchmark_file(self, template):
        parent, rel_path = self._get_relative_path(template)
        rel_path = rel_path.with_suffix(".log")
        return self.benchmark_dir.joinpath(parent, rel_path)

    def get_log_file(self, template):
        parent, rel_path = self._get_relative_path(template)
        rel_path = rel_path.with_suffix(".log")
        return self.log_dir.joinpath(parent, rel_path)

    @staticmethod
    def get_cfclone_run_type_args(wildcards):
        run_type = wildcards.run_type
        if run_type == "full":
            return ""
        elif run_type == "normal":
            return "--only-normal"
        else:
            clone = run_type.split("_")[-1]
            return "--use-clone {}".format(clone)

    def _get_ctdna_file(self, sample):
        return self.config["ctdna_files"][sample]

    def _get_relative_path(self, template):
        try:
            rel_path = template.relative_to(self.pipeline_dir)
            parent = "working"
        except ValueError:
            rel_path = template.relative_to(self.out_dir)
            parent = "output"
        return parent, rel_path
