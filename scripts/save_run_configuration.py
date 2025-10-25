import sys
import yaml


def save_run_config(config, config_yaml):
    with open(config_yaml, 'w', encoding='utf8') as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)


if __name__ == '__main__':
    sys.stderr = open(snakemake.log[0], "w")

    sys.stdout = sys.stderr

    save_run_config(snakemake.params["run_config"], snakemake.output["run_config"])
