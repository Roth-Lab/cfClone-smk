from snakemake.io import expand, directory
from snakemake.utils import min_version, validate

min_version("9.4")


conda: "envs/global.yaml"


validate(config, "schemas/config.schema.yaml")

from utils import ConfigManager

config = ConfigManager(config)


rule all:
    input:
        config.pipeline_files,
    localrule: True


rule run_cfclone:
    input:
        c=config.clone_cn_file,
        i=config.get_ctdna_file,
    output:
        config.fit_template,
    params:
        c=config.num_chains,
        r=config.num_rounds,
        v=config.num_chains_vi,
        rt=config.get_cfclone_run_type_args,
    benchmark:
        config.get_benchmark_file(config.fit_template)
    conda:
        "envs/cfclone.yaml"
    log:
        config.get_log_file(config.fit_template),
    threads: config.num_threads
    shell:
        "(cfclone fit "
        "-c {input.c} "
        "-i {input.i} "
        "-o {output} "
        "-t {threads} "
        "--num-chains {params.c} "
        "--num-chains-vi {params.v} "
        "--num-rounds {params.r} "
        "{params.rt})  >{log} 2>&1"


rule post_process_result:
    input:
        config.fit_template,
    output:
        config.post_process_template,
    params:
        lambda wildcards: wildcards.pp_result.replace("_", "-"),
    conda:
        "envs/cfclone.yaml"
    log:
        config.get_log_file(config.post_process_template),
    shell:
        "(cfclone write-{params} -i {input} -o {output}) >{log} 2>&1"


# rule all:
#     input:
#         expand(config_manager.prev_sized_nodes_tree_json, sample=config_manager.sample_list),
#         expand(config_manager.labelled_nodes_tree_json, sample=config_manager.sample_list),
#         expand(config_manager.prev_sized_nodes_tree_svg, sample=config_manager.sample_list),
#         expand(config_manager.labelled_nodes_tree_png, sample=config_manager.sample_list),
#         expand(config_manager.error_bars_tree_plot, sample=config_manager.sample_list),
#         expand(config_manager.rdr_baf_fit_plot, sample=config_manager.sample_list),
#         expand(config_manager.sample_rho_dominance_heatmap, sample=config_manager.sample_list),
#         config_manager.fits_plot,
#         config_manager.gathered_error_bars_tree_plot,
#         config_manager.clonal_prevalence_table,
#         config_manager.experiment_configuration,
#         config_manager.rho_dominance_plot,
#         config_manager.all_pigeons_summary_table,
#     localrule: True

# rule run_cfclone:
#     input:
#         c=config_manager.clone_cn_file,
#         i=config_manager.ctdna_file,
#     output:
#         config_manager.fit_file

# rule format_cn_profiles_from_bed:
#     input:
#         bed_file=config_manager.patient_cn_bed_file,
#     output:
#         formatted_cn_profiles=config_manager.patient_cn_profiles,
#     log:
#         config_manager.get_log_file("", "format_cn_profiles_from_bed.log"),
#     conda:
#         "envs/python.yaml"
#     group:
#         "pre-proc"
#     script:
#         "scripts/get_cfClone_cn_profiles_from_clone_bed.py"


# rule create_sample_cfClone_input_files:
#     input:
#         raw_ctDNA_file=config_manager.get_sample_raw_ctDNA_file,
#         patient_cn_profiles=config_manager.patient_cn_profiles,
#     output:
#         ctDNA_file=config_manager.sample_ctDNA_file,
#         cn_profile=config_manager.sample_cn_profile,
#     log:
#         config_manager.get_log_file("create_sample_cfClone_input_files", "{sample}.log"),
#     conda:
#         "envs/python.yaml"
#     group:
#         "sample-pre-proc"
#     script:
#         "scripts/format_cfClone_ctDNA_data.py"


# rule preprocess_clone_tree_create_colour_dict:
#     input:
#         clone_tree=config_manager.raw_clone_tree,
#         patient_cn_profiles=config_manager.patient_cn_profiles,
#     output:
#         processed_clone_tree=config_manager.labelled_newick,
#         clone_colour_dict_json=config_manager.colour_dict_json,
#     log:
#         config_manager.get_log_file("", "preprocess_clone_tree_create_colour_dict.log"),
#     conda:
#         "envs/plot.yaml"
#     group:
#         "pre-proc"
#     script:
#         "scripts/preprocess_input_tree.py"


# rule build_restricted_model_input:
#     input:
#         cn_profile=config_manager.sample_cn_profile,
#     output:
#         cn_profile=config_manager.restricted_sample_cn_profile,
#     params:
#         clone=lambda wildcards: wildcards.clone,
#     conda:
#         "envs/python.yaml"
#     log:
#         config_manager.get_log_file("build_restricted_model_input", "{sample}_{clone}.log"),
#     group:
#         "sample-clone-pre-proc"
#     script:
#         "scripts/build_restricted_model_input.py"


# rule run_cfclone:
#     input:
#         cn_profile=config_manager.sample_cn_profile,
#         ctDNA_file=config_manager.sample_ctDNA_file,
#     output:
#         results_dir=directory(config_manager.cfclone_sample_results_dir),
#         samples_file=config_manager.cfclone_results_file,
#         pigeons_summary=config_manager.raw_pigeons_summary,
#     threads: config_manager.cfclone_threads
#     conda:
#         "envs/cfclone.yaml"
#     params:
#         csv_path=config_manager.csv_path,
#         pigeons_summary=config_manager.raw_pigeons_summary_path,
#         subsampling=config_manager.subsampling,
#         n_rounds=config_manager.n_rounds,
#         symmetrization=config_manager.symmetrization,
#         baf_likelihood=config_manager.baf_likelihood,
#         rdr_likelihood=config_manager.rdr_likelihood,
#         keep_only_normal=config_manager.keep_only_normal,
#     log:
#         infer_log=config_manager.get_log_file("run_cfclone", "{sample}_infer.log"),
#         report_log=config_manager.get_log_file("run_cfclone", "{sample}_report.log"),
#     benchmark:
#         config_manager.get_benchmark_file("run_cfclone", "{sample}.txt")
#     group:
#         "run-cfClone"
#     shell:
#         """
#         mkdir -p {output.results_dir}

#         cd {output.results_dir}

#         (cfclone infer \
#         --cfdna-tsv {input.ctDNA_file} \
#         --clone-cn-tsv {input.cn_profile} \
#         --subsampling {params.subsampling} \
#         --n-rounds {params.n_rounds} \
#         --symmetrization {params.symmetrization} \
#         --baf-likelihood {params.baf_likelihood} \
#         --rdr-likelihood {params.rdr_likelihood} \
#         --n-threads {threads} {params.keep_only_normal}) >{log.infer_log} 2>&1

#         (cfclone run-report --save-csv) >{log.report_log} 2>&1

#         cp {params.csv_path} {output.samples_file}

#         cp {params.pigeons_summary} {output.pigeons_summary}
#         """


# use rule run_cfclone as run_cfclone_restricted_model with:
#     input:
#         cn_profile=config_manager.restricted_sample_cn_profile,
#         ctDNA_file=config_manager.sample_ctDNA_file,
#     output:
#         results_dir=directory(config_manager.restricted_cfclone_sample_results_dir),
#         samples_file=config_manager.restricted_cfclone_results_file,
#         pigeons_summary=config_manager.raw_restricted_pigeons_summary,
#     log:
#         infer_log=config_manager.get_log_file("run_cfclone", "{sample}_{clone}_infer.log"),
#         report_log=config_manager.get_log_file("run_cfclone", "{sample}_{clone}_report.log"),
#     benchmark:
#         config_manager.get_benchmark_file("run_cfclone", "{sample}_{clone}.txt")
#     group:
#         "run-restricted-cfClone"


# rule process_pigeons_summary:
#     input:
#         raw_pigeons_summary=config_manager.raw_pigeons_summary,
#     output:
#         pigeons_summary=config_manager.pigeons_summary,
#     conda:
#         "envs/python.yaml"
#     log:
#         config_manager.get_log_file("process_pigeons_summary", "{sample}.log"),
#     group:
#         "run-cfClone"
#     params:
#         model_used="full_model",
#         sample_id=lambda wildcards: wildcards.sample,
#     script:
#         "scripts/process_pigeons_summary.py"


# use rule process_pigeons_summary as process_restricted_model_pigeons_summary with:
#     input:
#         raw_pigeons_summary=config_manager.raw_restricted_pigeons_summary,
#     output:
#         pigeons_summary=config_manager.restricted_pigeons_summary,
#     log:
#         config_manager.get_log_file("process_pigeons_summary", "{sample}_{clone}.log"),
#     group:
#         "run-restricted-cfClone"
#     params:
#         model_used=config_manager.get_restricted_model_name,
#         sample_id=lambda wildcards: wildcards.sample,


# rule gather_sample_pigeons_summaries:
#     input:
#         expand(
#             config_manager.restricted_pigeons_summary,
#             sample=lambda wildcards: wildcards.sample,
#             clone=config_manager.clones,
#         ),
#         config_manager.pigeons_summary,
#     output:
#         out_table=config_manager.sample_pigeons_summary,
#     log:
#         config_manager.get_log_file("gather_sample_pigeons_summaries", "{sample}.log"),
#     script:
#         "scripts/gather_tables.py"


# use rule gather_sample_pigeons_summaries as gather_all_pigeons_summaries with:
#     input:
#         expand(config_manager.sample_pigeons_summary, sample=config_manager.sample_list),
#     output:
#         out_table=config_manager.all_pigeons_summary_table,
#     log:
#         config_manager.get_log_file("", "gather_all_pigeons_summaries.log"),


# rule create_clonal_rho_table_from_cfClone_result:
#     input:
#         cfClone_samples=config_manager.cfclone_results_file,
#     output:
#         cfClone_clonal_rho_table=config_manager.cfclone_samples_clonal_rho,
#     conda:
#         "envs/python.yaml"
#     log:
#         config_manager.get_log_file("create_clonal_rho_table_from_cfClone_result", "{sample}.log"),
#     group:
#         "sample-post-proc"
#     benchmark:
#         config_manager.get_benchmark_file("create_clonal_rho_table_from_cfClone_result", "{sample}.txt")
#     script:
#         "scripts/cfClone_samples_convert_to_long.py"


# rule create_mu_p_tables_from_cfClone_result:
#     input:
#         cfClone_samples=config_manager.cfclone_results_file,
#     output:
#         mu_table=config_manager.cfClone_mu_table,
#         p_table=config_manager.cfClone_p_table,
#     conda:
#         "envs/plot.yaml"
#     log:
#         config_manager.get_log_file("create_mu_p_tables_from_cfClone_result", "{sample}.log"),
#     group:
#         "sample-post-proc"
#     benchmark:
#         config_manager.get_benchmark_file("create_mu_p_tables_from_cfClone_result", "{sample}.txt")
#     script:
#         "scripts/samples_convert_to_long_mu_p.py"


# rule create_clonal_prev_summary_table_and_tree_plots:
#     input:
#         clonal_rho_table=config_manager.cfclone_samples_clonal_rho,
#         processed_clone_tree=config_manager.labelled_newick,
#         clone_colour_dict_json=config_manager.colour_dict_json,
#     output:
#         rho_summary_table=config_manager.sample_clonal_rho_summary,
#         prev_sized_nodes_tree_json=config_manager.prev_sized_nodes_tree_json,
#         labelled_nodes_tree_json=config_manager.labelled_nodes_tree_json,
#         prev_sized_nodes_tree_svg=config_manager.prev_sized_nodes_tree_svg,
#         labelled_nodes_tree_png=config_manager.labelled_nodes_tree_png,
#     log:
#         config_manager.get_log_file("create_clonal_prev_summary_table_and_tree_plots", "{sample}.log"),
#     params:
#         sample_id=lambda wildcards: wildcards.sample,
#     conda:
#         "envs/plot.yaml"
#     group:
#         "sample-post-proc"
#     benchmark:
#         config_manager.get_benchmark_file("create_clonal_prev_summary_table_and_tree_plots", "{sample}.txt")
#     script:
#         "scripts/create_rho_summary.py"


# rule create_sample_error_bars_tree_plot:
#     input:
#         rho_summary_table=config_manager.sample_clonal_rho_summary,
#         prev_sized_nodes_tree_json=config_manager.prev_sized_nodes_tree_json,
#     output:
#         error_bars_tree_plot=config_manager.error_bars_tree_plot,
#     log:
#         config_manager.get_log_file("create_error_bars_tree_plot", "{sample}.log"),
#     conda:
#         "envs/plot.yaml"
#     group:
#         "sample-plotting"
#     script:
#         "scripts/error_bars_tree_plot.py"
# rule create_sample_fit_plot:
#     input:
#         mu_table=config_manager.cfClone_mu_table,
#         p_table=config_manager.cfClone_p_table,
#         ctDNA_file=config_manager.sample_ctDNA_file,
#     output:
#         fit_plot=config_manager.rdr_baf_fit_plot,
#     log:
#         config_manager.get_log_file("create_sample_fit_plot", "{sample}.log"),
#     params:
#         sample_id=lambda wildcards: wildcards.sample,
#     conda:
#         "envs/plot.yaml"
#     group:
#         "sample-plotting"
#     script:
#         "scripts/plot_cfClone_fits.py"
# rule create_sample_rho_dominance_heatmap:
#     input:
#         clonal_rho_table=config_manager.cfclone_samples_clonal_rho,
#         processed_clone_tree=config_manager.labelled_newick,
#     output:
#         rho_heatmap=config_manager.sample_rho_dominance_heatmap,
#     log:
#         config_manager.get_log_file("create_sample_rho_dominance_tree_heatmap", "{sample}.log"),
#     params:
#         sample_id=lambda wildcards: wildcards.sample,
#     conda:
#         "envs/plot.yaml"
#     group:
#         "sample-plotting"
#     script:
#         "scripts/rho_dominance_heatmap.py"
# rule create_gathered_fits_plot:
#     input:
#         mu_tables=expand(config_manager.cfClone_mu_table, sample=config_manager.sample_list),
#         p_tables=expand(config_manager.cfClone_p_table, sample=config_manager.sample_list),
#         ctDNA_files=expand(config_manager.sample_ctDNA_file, sample=config_manager.sample_list),
#     output:
#         fit_plot=config_manager.fits_plot,
#     params:
#         sample_ids=config_manager.sample_list,
#     log:
#         config_manager.get_log_file("", "create_gathered_fits_plot.log"),
#     conda:
#         "envs/plot.yaml"
#     group:
#         "gathered-plotting"
#     script:
#         "scripts/plot_all_sample_fits.py"
# rule create_gathered_error_bars_tree_plot:
#     input:
#         rho_tables=expand(config_manager.sample_clonal_rho_summary, sample=config_manager.sample_list),
#         tree_jsons=expand(config_manager.prev_sized_nodes_tree_json, sample=config_manager.sample_list),
#     output:
#         errorbar_tree_plot=config_manager.gathered_error_bars_tree_plot,
#     log:
#         config_manager.get_log_file("", "create_gathered_error_bars_tree_plot.log"),
#     conda:
#         "envs/plot.yaml"
#     group:
#         "gathered-plotting"
#     script:
#         "scripts/plot_all_sample_error_bar_trees.py"
# rule create_gathered_clonal_prev_stats_table:
#     input:
#         rho_tables=expand(config_manager.sample_clonal_rho_summary, sample=config_manager.sample_list),
#     output:
#         clonal_prevalence_table=config_manager.clonal_prevalence_table,
#     log:
#         config_manager.get_log_file("", "create_gathered_clonal_prev_stats_table.log"),
#     conda:
#         "envs/python.yaml"
#     group:
#         "gathered-plotting"
#     script:
#         "scripts/gathered_clonal_prev_table.py"
# rule create_gathered_rho_dominance_heatmap:
#     input:
#         rho_tables=expand(config_manager.cfclone_samples_clonal_rho, sample=config_manager.sample_list),
#         tree_newick=config_manager.labelled_newick,
#     output:
#         rho_dominance_plots=config_manager.rho_dominance_plot,
#     log:
#         config_manager.get_log_file("", "create_gathered_rho_dominance_heatmap.log"),
#     conda:
#         "envs/plot.yaml"
#     group:
#         "gathered-plotting"
#     params:
#         sample_ids=config_manager.sample_list,
#     script:
#         "scripts/plot_all_sample_rho_dominance_heatmaps.py"
# rule save_run_configuration:
#     output:
#         run_config=config_manager.experiment_configuration,
#     log:
#         config_manager.get_log_file("", "save_run_configuration.log"),
#     conda:
#         "envs/python.yaml"
#     params:
#         run_config=config_manager.config,
#     group:
#         "pre-proc"
#     script:
#         "scripts/save_run_configuration.py"
