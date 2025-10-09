import pandas as pd
from skbio.tree import TreeNode
import networkx as nx
from seaborn import color_palette
import colorcet
from json_serialisation import serialise_dict_to_json
import sys


def preprocess_tree_to_extant_tips(df, tree):
    index_set = set(df["clone"].unique())
    tip_set = set(tip.name for tip in tree.tips())
    tip_intersect = tip_set.intersection(index_set)
    return tree.shear(tip_intersect, prune=True, strict=False, inplace=False)


def assign_ancestral_node_names(df, treefile):
    tree = TreeNode.read(treefile, convert_underscores=False)
    tree = preprocess_tree_to_extant_tips(df, tree)
    template = "ancestral_{}"
    tree.assign_ids()

    ancestral_id_map = {}

    inner_node_idx = 0

    node_names = []

    for node in tree.preorder():
        if node.is_tip():
            continue
        ancestral_id_map[node.id] = inner_node_idx
        inner_node_idx += 1

    for node in tree.postorder():

        if node.name is None:
            node.name = template.format(ancestral_id_map[node.id])

        node_names.append(node.name)

    return tree, sorted(node_names)


def create_clone_colour_dict(node_names):
    num_colours = len(node_names)
    if num_colours > 10:
        colour_list = color_palette(colorcet.glasbey, num_colours).as_hex()
    else:
        colour_list = color_palette('colorblind', num_colours).as_hex()
    gt_colour_dict = dict(zip(node_names, colour_list))
    return gt_colour_dict


def main(in_tree_newick, cn_profiles_df, out_tree_newick, out_colour_dict_json):
    df = pd.read_table(cn_profiles_df)
    labelled_tree, node_names = assign_ancestral_node_names(df, in_tree_newick)
    clone_colour_dict = create_clone_colour_dict(node_names)
    labelled_tree.write(out_tree_newick, format="newick")
    serialise_dict_to_json(clone_colour_dict, out_colour_dict_json)


if __name__ == '__main__':
    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    main(
        snakemake.input["clone_tree"],
        snakemake.input["patient_cn_profiles"],
        snakemake.output["processed_clone_tree"],
        snakemake.output["clone_colour_dict_json"],
    )
