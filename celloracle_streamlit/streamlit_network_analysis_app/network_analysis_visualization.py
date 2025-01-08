# Streamlit practice 2021 0223

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns




def plot_score_comparison(merged_score, value, cluster1, cluster2):

    import plotly.express as px
    res = merged_score[merged_score.cluster.isin([cluster1, cluster2])][[value, "cluster"]]
    res = res.reset_index(drop=False)
    piv = pd.pivot_table(res, values=value, columns="cluster", index="index")
    piv = piv.reset_index(drop=False)
    piv = piv.fillna(0)

    fig = px.scatter(piv, x=cluster1, y=cluster2,
                     hover_data=['index'], template="plotly_white")

    return fig



def plot_carto_terms(merged_score, palette, gene):

    fig = plt.figure(figsize=[3, 6])

    print(gene)
    tt = pd.get_dummies(merged_score[["cluster", "role"]],columns=["role"])
    tt = tt.loc[gene].set_index("cluster")
    tt.columns = [i.replace("role_", "") for i in tt.columns]
    order = ["Ultra peripheral", "Peripheral", "Connector","Kinless","Provincical Hub","Connector Hub", "Kinless Hub"]
    tt = tt.reindex(index=palette.index.values, columns=order).astype("float").fillna(0)

    sns.heatmap(data=tt, cmap="Blues", cbar=False)

    return fig


def plot_carto(merged_score, gene, cluster):
    data = merged_score[merged_score.cluster == cluster]
    data = data[["connectivity", "participation"]]
    data = data.dropna(axis=0)

    fig = plt.figure()

    _plot_base_cartography(z_min=-1, z_max=data.connectivity.max())

    plt.scatter(data.participation.values,
                data.connectivity.values,
                marker='o', edgecolor="lightgray",
                c="none", alpha=1)

    if gene in data.index:
        x = data.loc[gene, "participation"]
        y = data.loc[gene, "connectivity"]
        args_annot = {"size": 10, "color": "black"}
        arrow_dict = {"width": 0.5, "headwidth": 0.5, "headlength": 1, "color": "black"}
        plt.annotate(gene, xy=(x, y), xytext=(x+0.05, y+0.05),
                     arrowprops=arrow_dict, **args_annot)

    plt.xlabel("Participation coefficient (P)")
    plt.ylabel("Whithin-module\ndegree (z)")

    return fig


def _plot_base_cartography(z_min, z_max):
    args_line = {"linestyle": "dashed",
               "alpha": 0.5,
               "c": "gray"}
    plt.plot([-0.05, 1.05], [2.5, 2.5], **args_line)
    plt.plot([0.05, 0.05], [z_min, 2.5], **args_line)
    plt.plot([0.62, 0.62], [z_min, 2.5], **args_line)
    plt.plot([0.8, 0.8], [z_min, 2.5], **args_line)
    plt.plot([0.3, 0.3], [2.5, z_max + 0.5], **args_line)
    plt.plot([0.75, 0.75], [2.5, z_max + 0.5], **args_line)
    plt.xlabel("Participation coefficient (P)")
    plt.ylabel("Whithin-module degree (z)")
    plt.xlim([-0.1, 1.1])

def plot_scores_as_rank(merged_score, cluster, value, n_gene=50):


    res = merged_score[merged_score.cluster == cluster]
    res = res[value].sort_values(ascending=False)
    res = res[:n_gene]

    fig = plt.figure(figsize=[3, 6])

    plt.scatter(res.values, range(len(res)))
    plt.yticks(range(len(res)), res.index.values)#, rotation=90)
    #plt.xlabel(value)
    plt.title(f" {value} \n top {n_gene} in {cluster}")
    plt.ticklabel_format(style='sci',
                         axis='x',
                         scilimits=(0,0))

    plt.gca().invert_yaxis()


    #plt.subplots_adjust(left=0.5, right=0.99)

    return fig

def add_label_cartography(data):

    # categorize nodes
    R1_nodes = data[(data.participation <= 0.05) & (data.connectivity < 2.5)].index
    R2_nodes = data[(data.participation > 0.05) & (data.participation <= 0.62) & (data.connectivity < 2.5)].index
    R3_nodes = data[(data.participation > 0.62) & (data.participation <= 0.80) & (data.connectivity < 2.5)].index
    R4_nodes = data[(data.participation > 0.80) & (data.connectivity < 2.5)].index
    R5_nodes = data[(data.participation <= 0.30) & (data.connectivity >= 2.5)].index
    R6_nodes = data[(data.participation > 0.30) & (data.participation <= 0.75) & (data.connectivity >= 2.5)].index
    R7_nodes = data[(data.participation > 0.75) & (data.connectivity >= 2.5)].index

    node_kind = ["R1: Ultra peripheral", "R2: Peripheral", "R3: Non-hub",
                 "R4: Non-hub kinless", "R5: Provincial",
                 "R6: Connector hubs", "R7: Kinless hubs"]

    data["label_cartography"] = "0"

    for i in range(1, 8):
        data.loc[eval(f"R{i}_nodes"), "label_cartography"] = node_kind[(i-1)]

    return data


def plot_score_per_cluster(merged_score, palette, gene, figsize=[6, 5], plt_show=False):
    """
    Plot network score for a specific gene.
    This function can be used to compare network score of a specific gene between clusters
    and get insight about the dynamics of the gene.

    Args:
        links (Links object): See network_analisis.Links class for detail.
        gene (srt): Gene name.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """
    fig = plt.figure(figsize=figsize)
    print(gene)
    res = merged_score[merged_score.index==gene]
    res = res.rename(
        columns={"degree_centrality_all": "degree\ncentrality",
                 "betweenness_centrality": "betweenness\ncentrality",
                 "eigenvector_centrality": "eigenvector\ncentrality"})
    # make plots
    values = [ "degree\ncentrality",  "betweenness\ncentrality",
              "eigenvector\ncentrality"]
    for i, value in zip([1, 2, 3], values):
        plt.subplot(1, 3,  i)
        ax = sns.stripplot(data=res, y="cluster", x=value,
                      size=10, orient="h",linewidth=1, edgecolor="w",
                      order=palette.index.values,
                      palette=palette.palette.values)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.yaxis.grid(True)
        ax.tick_params(bottom=False,
                        left=False,
                        right=False,
                        top=False)
        if i > 1:
            plt.ylabel(None)
            ax.tick_params(labelleft=False)

    return fig
