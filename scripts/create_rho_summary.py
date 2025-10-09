import pandas as pd
from skbio.tree import TreeNode
import xarray as xr
import arviz as az
import networkx as nx
from seaborn import color_palette
import colorcet
import math
from json_serialisation import serialise_networkx_tree_to_json, read_json_into_dict
import sys


def preprocess_tree_to_extant_tips(df, tree):
    index_set = set(df["clone_id"].unique())
    tip_set = set(tip.name for tip in tree.tips())
    tip_intersect = tip_set.intersection(index_set)
    return tree.shear(tip_intersect, prune=True, strict=False, inplace=False)


def draw_tree(tree, filename, transparent_bg=True):
    draw_graph = nx.nx_agraph.to_agraph(tree)

    draw_graph.graph_attr['dpi'] = '300'
    if transparent_bg:
        draw_graph.graph_attr['bgcolor'] = "transparent"
    draw_graph.graph_attr['fontname'] = 'Nimbus Sans'
    draw_graph.node_attr['fontname'] = 'Nimbus Sans'
    draw_graph.edge_attr['arrowhead'] = 'normal'

    draw_graph.draw(filename, prog="dot")


def scikit_tree_to_networkx_add_prev_and_hdi(tree, df_summary, colour_dict, prev_size_nodes=True, ancestors_grey=True):
    nx_graph = nx.DiGraph()

    index_template = 'rho[{}]'

    tip_label = "Clone {}\nClonal Prev: {}\nHDI lower: {} | HDI upper: {}"
    inner_label = "{}\nClonal Prev: {}\nHDI lower: {} | HDI upper: {}"

    base_size = 1
    base_radius = base_size / 2
    base_area = math.pi * (base_radius**2)

    grey_to_use = "#696969"

    for node in tree.traverse():

        node_name = node.name
        summary_row = df_summary.loc[index_template.format(node_name)]

        for child in node.children:
            child_name = child.name
            nx_graph.add_edge(node_name, child_name, length=child.length)

        hdi_low = summary_row["hdi_lower"]
        hdi_upper = summary_row["hdi_upper"]
        mean = summary_row["mean"]

        if node.is_tip():
            label = tip_label.format(node_name, round(mean, 3), round(hdi_low, 3), round(hdi_upper, 3))
            node_colour = colour_dict[node_name]
        else:
            label = inner_label.format(node_name, round(mean, 3), round(hdi_low, 3), round(hdi_upper, 3))
            node_colour = colour_dict[node_name]
            if ancestors_grey:
                node_colour = grey_to_use

        if prev_size_nodes:
            if node_colour == grey_to_use:
                fillcolor = node_colour + '33'
            else:
                fillcolor = node_colour + '33'
            area_covered = base_area * mean
            node_radius = math.sqrt(area_covered / math.pi)
            size_mod = node_radius * 2
            node_attr_dict = {
                "label": "",
                "shape": "circle",
                "width": size_mod,
                "fixedsize": True,
                "style": "filled",
                "penwidth": 0.75,
                "color": node_colour,
                "fillcolor": fillcolor,
            }
        else:
            if node_colour == grey_to_use:
                fillcolor = "white"
            else:
                fillcolor = node_colour + '1A'
            node_attr_dict = {
                "label": label,
                "color": node_colour,
                "fillcolor": fillcolor,
                "style": "filled",
                "penwidth": 3,
            }

        # hdi_width = summary_row["hdi_width_normed"]
        # hdi_alpha_pct = hdi_width
        # hdi_alpha = '{:02x}'.format(round(255 * hdi_alpha_pct))

        nx_graph.add_node(node_name, **node_attr_dict, rho_stats=summary_row)

    return nx_graph


def make_clonal_prev_tree_figure(
    clonal_rho_table,
    prev_sized_tree_image,
    treefile,
    out_summary_df,
    prev_sized_tree_json,
    sample_id,
    labelled_nodes_tree_json,
    labelled_nodes_tree_image,
    colour_dict_json,
):

    clones_to_remove = {"normal"}

    df = pd.read_table(clonal_rho_table, usecols=["iteration", "chain", "clone_id", "rho"])

    df = remove_clones_normalize_rho(clones_to_remove, df)

    df, tree = build_clone_df_and_tree(df, treefile)

    colour_dict = read_json_into_dict(colour_dict_json)

    df_summary = build_arviz_rho_summary_df(df)

    define_hdi_upper_and_lower_cols(df_summary)

    # compute_hdi_width_cols(df_summary)

    nx_tree = scikit_tree_to_networkx_add_prev_and_hdi(
        tree,
        df_summary,
        colour_dict,
        prev_size_nodes=True,
        ancestors_grey=True,
    )
    nx_tree.graph["sample_id"] = sample_id
    serialise_networkx_tree_to_json(nx_tree, prev_sized_tree_json)

    draw_tree(nx_tree, prev_sized_tree_image, transparent_bg=True)

    nx_tree = scikit_tree_to_networkx_add_prev_and_hdi(
        tree,
        df_summary,
        colour_dict,
        prev_size_nodes=False,
        ancestors_grey=True,
    )
    nx_tree.graph["sample_id"] = sample_id
    serialise_networkx_tree_to_json(nx_tree, labelled_nodes_tree_json)

    draw_tree(nx_tree, labelled_nodes_tree_image, transparent_bg=False)

    df_summary = finalise_rho_summary_df(colour_dict, df_summary, sample_id, tree)

    df_summary.to_csv(out_summary_df, sep="\t", index=False)


def finalise_rho_summary_df(colour_dict, df_summary, sample_id, tree):
    df_summary["sample_id"] = sample_id
    df_summary = df_summary.reset_index(names=["parameters"])
    df_summary[["parameter_name", "clone_id"]] = df_summary["parameters"].str.split('[', expand=True)
    df_summary["clone_id"] = df_summary["clone_id"].str.removesuffix(']')
    node_order_dict = {v.name: k for k, v in enumerate(tree.postorder())}
    df_summary = df_summary.set_index("clone_id")
    df_summary["clone_tree_order_idx"] = node_order_dict
    df_summary["clone_colour"] = colour_dict
    df_summary = df_summary.reset_index()
    return df_summary


def compute_hdi_width_cols(df_summary):
    df_summary["hdi_width"] = df_summary["hdi_upper"] - df_summary["hdi_lower"]
    width_sum = df_summary["hdi_width"].max()
    df_summary["hdi_width_normed"] = df_summary["hdi_width"] / width_sum


def define_hdi_upper_and_lower_cols(df_summary):
    hdi_cols = {col: float(col[4:-1]) for col in df_summary.columns if col.startswith("hdi")}
    hdi_col_names = list(hdi_cols.keys())
    if hdi_cols[hdi_col_names[0]] < hdi_cols[hdi_col_names[1]]:
        hdi_col_name_map = {hdi_col_names[0]: "hdi_lower", hdi_col_names[1]: "hdi_upper"}
    else:
        hdi_col_name_map = {hdi_col_names[1]: "hdi_lower", hdi_col_names[0]: "hdi_upper"}
    for hdi_col, new_name in hdi_col_name_map.items():
        df_summary[new_name] = df_summary[hdi_col]


def build_arviz_rho_summary_df(df):
    df = df.rename(columns={"iteration": "draw"})
    df = df.set_index(["chain", "draw", "clone_id"])
    xdata = xr.Dataset.from_dataframe(df)
    az_dataset = az.InferenceData(posterior=xdata)
    df_summary = az.summary(
        az_dataset,
        var_names=["rho"],
        kind="stats",
        hdi_prob=0.95,
        round_to="none",
    )
    return df_summary


def remove_clones_normalize_rho(clones_to_remove, df):
    df = df.loc[~df["clone_id"].isin(clones_to_remove)]
    grouped = df.groupby(["iteration", "chain"])
    normed_grps = []
    for name, group in grouped:
        grp_sum = group["rho"].sum()
        group["rho"] /= grp_sum
        normed_grps.append(group)
    df = pd.concat(normed_grps)
    return df


def build_clone_df_and_tree(df, treefile):
    tree = TreeNode.read(treefile, convert_underscores=False)
    tree = preprocess_tree_to_extant_tips(df, tree)
    df = df.set_index(["iteration", "chain"])
    clonal_grouped = df.groupby("clone_id")

    grp_dict = dict()
    for node in tree.postorder():
        if node.is_tip():
            grp_dict[node.name] = clonal_grouped.get_group(node.name)
            continue

        node_grp = grp_dict[node.children[0].name].copy()
        node_grp["clone_id"] = node.name

        for child in node.children[1:]:
            child_name = child.name
            child_grp = grp_dict[child_name]
            node_grp["rho"] += child_grp["rho"]

        grp_dict[node.name] = node_grp
    clone_df = pd.concat(grp_dict.values()).reset_index()
    return clone_df, tree


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    sys.stdout = sys.stderr

    make_clonal_prev_tree_figure(
        snakemake.input["clonal_rho_table"],
        snakemake.output["prev_sized_nodes_tree_svg"],
        snakemake.input["processed_clone_tree"],
        snakemake.output["rho_summary_table"],
        snakemake.output["prev_sized_nodes_tree_json"],
        snakemake.params["sample_id"],
        snakemake.output["labelled_nodes_tree_json"],
        snakemake.output["labelled_nodes_tree_png"],
        snakemake.input["clone_colour_dict_json"],
    )
