# Streamlit practice 2021 0223

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

import celloracle as co
from celloracle.network_analysis.cartography import plot_cartography_kde


def network_analysis_01(path_links, embedding_key,
    default_cluster_0, default_cluster_1, cluster_column_name,
    default_gene,
    network_score_kinds, n_genes,
    description_0, title, path_adata=None):


    ## Define functions and hyperoarameters
    @st.cache(allow_output_mutation=True)
    def load_anndata(path_adata):
        #path = "test3.h5ad"
        adata = sc.read_h5ad(path_adata)

        return adata

    @st.cache(allow_output_mutation=True)
    def load_links(path_links):
        #path = "test3.h5ad"
        links = co.load_hdf5(path_links)

        return links


    def plot_embeddings(x, args={}):
        fig = sc.pl.embedding(adata=adata, basis=embedding_key, color=x, s=50, return_fig=True, **args)
        return fig

    def plot_scores_as_rank(cluster, value, n_gene=50):


        res = links.merged_score[links.merged_score.cluster == cluster]
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


    def plot_score_per_cluster(gene):
        fig = plt.figure(figsize=[6,5])
        links.plot_score_per_cluster(goi=gene, plt_show=False)
        return fig


    def plot_score_comparison(value, cluster1, cluster2, percentile):
        fig = plt.figure(figsize=[5,5])
        links.plot_score_comparison_2D(value=value,
                                       cluster1=cluster1,
                                       cluster2=cluster2,
                                       percentile=percentile)
        return fig

    def plot_carto(gene, cluster):
        data = links.merged_score[links.merged_score.cluster == cluster]
        data = data[["connectivity", "participation"]]
        fig = plt.figure()
        plot_cartography_kde(data,
                             gois=gene,
                             scatter=True,
                             kde=False,
                             args_kde={"n_levels": 105},
                             args_line={"c":"gray"})

        return fig

    def plot_carto_terms(gene):
        fig = plt.figure(figsize=[3, 6])
        links.plot_cartography_term(goi=gene, plt_show=False)
        #plt.yticks(links.cluster, fontsize=5)
        return fig


    ### APP starts here

    ## Load data
    if path_adata is not None:
        adata = load_anndata(path_adata=path_adata)
    links = load_links(path_links=path_links)

    genes_in_links = links.merged_score.index.unique()

    ## Side bar

    st.sidebar.write(f"# {title}")

    st.sidebar.write("Please select data below.")

    # Select cluster
    st.sidebar.write("## Cluster")
    clusters = links.cluster
    cluster = st.sidebar.selectbox("Select cluster to show detailed data",
                                   clusters, index=clusters.index(default_cluster_0))


    st.sidebar.write("## Gene")
    # Select gene
    genes = sorted(adata.var.index.values)
    gene = st.sidebar.selectbox("Select gene to analyze",
                                genes,
                                index=genes.index(default_gene))



    st.sidebar.write("## Network score")
    # Select metric
    value = st.sidebar.selectbox("Select network score for comparison",
                                 network_score_kinds,
                                 index=network_score_kinds.index("degree_centrality_all"))


    st.sidebar.write("## Score comparison between two GRNs")
    # Select two clusters to compare
    cluster1 = st.sidebar.selectbox("Select cluster1 for comparison",
                                   clusters, index=clusters.index(default_cluster_0))

    cluster2 = st.sidebar.selectbox("Select cluster2 for comparison",
                                   clusters, index=clusters.index(default_cluster_1))

    ## Main

    st.write(f"# {title}")
    st.write("## About")

    st.write(description_0)

    if path_adata is not None:
        st.write("# scRNA-seq data")

        #col1, col2 = st.beta_columns(2)
        st.write("### Clustering")
        st.pyplot(plot_embeddings(cluster_column_name, args={"legend_loc": "on data"}))

        #col2.write("### Gene expression")
        #col2.pyplot(plot_embeddings(gene, args={"use_raw":False, "cmap": "viridis"}))



    st.write(f"# Network scores")
    st.write(f"## 1. Top {n_genes} genes in {cluster} GRN centrality scores")

    col1, col2, col3 = st.beta_columns(3)
    col1.pyplot(plot_scores_as_rank(cluster=cluster,
                                    value="degree_centrality_all",
                                    n_gene=n_genes))
    col2.pyplot(plot_scores_as_rank(cluster=cluster,
                                    value="betweenness_centrality",
                                    n_gene=n_genes))
    col3.pyplot(plot_scores_as_rank(cluster=cluster,
                                    value="eigenvector_centrality",
                                    n_gene=n_genes))


    st.write(f"## 2. All scores in {cluster} GRN")

    df = links.merged_score[links.merged_score.cluster == cluster]
    st.dataframe(df)

    st.write(f"## 3. Network score comparison between GRNs")
    st.write(f"### GRN 1: {cluster1},     GRN 2: {cluster2}")
    st.write(f"### Network score : {value}")
    st.plotly_chart(
        links.plot_score_comparison_2D(
            value="degree_centrality_all",
            cluster1=cluster1, cluster2=cluster2,
            interactive=True))


    st.write("# Network dynamics")
    st.write(f"### Selected gene: {gene}")
    if gene in genes_in_links:
        pass
    else:
        st.write(f"### Caution! {gene} does not have any connection in the filtered GRNs.")
    st.write(f"## 1. Score dynamics")
    st.pyplot(plot_score_per_cluster(gene))

    st.write("## 2. Gene cartography analysis")
    col1, col2 = st.beta_columns([1, 2])
    col1.write("### 2.1 Cartography summary")
    if gene in genes_in_links:
        col1.pyplot(plot_carto_terms(gene=gene))
    else:
        col1.write(f"{gene} does not have any connection in the filtered GRNs.")
    col2.write(f"### 2.2 Cartography of {gene} in {cluster} GRN")
    col2.pyplot(plot_carto(gene=gene, cluster=cluster))
