from snakemake.utils import min_version, validate

min_version("9.4")

conda: "envs/global.yaml"

validate(config,"schemas/config.schema.yaml")

include: "utils.py"

config = ConfigManager(config)

ruleorder: merge_evidence > post_process_full_result


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


rule post_process_full_result:
    input:
        lambda wildcards: str(config.fit_template).format(sample=wildcards.sample,run_type="full"),
    output:
        config.post_process_template,
    params:
        lambda wildcards: wildcards.pp_result.replace("_","-"),
    conda:
        "envs/cfclone.yaml"
    log:
        config.get_log_file(config.post_process_template),
    shell:
        "(cfclone write-{params} -i {input} -o {output}) >{log} 2>&1"


rule write_ancestral_prevlances:
    input:
        i=lambda wildcards: str(config.fit_template).format(sample=wildcards.sample,run_type="full"),
        t=config.clone_tree_file,
    output:
        o=config.ancestral_prevalence_template,
        t=config.clone_prevalence_tree_template,
    conda:
        "envs/cfclone.yaml"
    log:
        config.get_log_file(config.ancestral_prevalence_template),
    shell:
        "(cfclone write-ancestral-prevalences -c {input.t} -i {input.i} -o {output.o} -t {output.t}) >{log} 2>&1"


rule write_evidence:
    input:
        config.fit_template,
    output:
        temp(config.run_type_evidence_template),
    conda:
        "envs/cfclone.yaml"
    log:
        config.get_log_file(config.run_type_evidence_template),
    shell:
        "(echo -n {wildcards.sample},{wildcards.run_type}, > {output}; "
        "cfclone print-model-evidence -i{input} >> {output}) 2>{log}"


rule merge_evidence:
    input:
        lambda wildcards: [
            str(config.run_type_evidence_template).format(run_type=r,sample=wildcards.sample)
            for r in config.run_types
        ],
    output:
        config.evidence_template,
    params:
        script=workflow.source_path("scripts/merge_evidence.py"),
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.evidence_template),
    shell:
        "(python {params.script} -i {input} -o {output}) >{log} 2>&1"


rule plot_clone_prevalences:
    input:
        d=config.ancestral_prevalence_template,
        t=config.clone_prevalence_tree_template,
    output:
        config.clone_prevalences_plot,
    params:
        script=workflow.source_path("scripts/plot_clone_prevalences.py"),
    conda:
        "envs/plot.yaml"
    log:
        config.get_log_file(config.clone_prevalence_tree_template),
    shell:
        "(python {params.script} -d {input.d} -t {input.t} -o {output}) >{log} 2>&1"


rule plot_fit:
    input:
        lambda wildcards: str(config.post_process_template).format(
            pp_result="parameter_summaries",run_type="full",sample=wildcards.sample
        ),
    output:
        config.fit_plot_template,
    params:
        script=workflow.source_path("scripts/plot_fit.py"),
    conda:
        "envs/plot.yaml"
    log:
        config.get_log_file(config.fit_plot_template),
    shell:
        "(python {params.script} -i {input} -o {output}) >{log} 2>&1"


rule plot_pairwsie_ranks:
    input:
        i=lambda wildcards: str(config.post_process_template).format(
            pp_result="pairwise_ranks",run_type="full",sample=wildcards.sample
        ),
        t=config.clone_tree_file,
    output:
        config.pairwise_ranks_plot,
    params:
        script=workflow.source_path("scripts/plot_pairwise_ranks.py"),
    conda:
        "envs/plot.yaml"
    log:
        config.get_log_file(config.pairwise_ranks_plot),
    shell:
        "(python {params.script} -i {input.i} -t {input.t} -o {output}) >{log} 2>&1"


rule save_run_configuration:
    output:
        run_config=config.experiment_configuration,
    log:
        config.log_dir.joinpath("save_run_configuration.log"),
    conda:
        "envs/python.yaml"
    params:
        run_config=config.config,
    script:
        "scripts/save_run_configuration.py"
