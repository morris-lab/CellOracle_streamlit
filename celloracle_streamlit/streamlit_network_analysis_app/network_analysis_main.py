# Streamlit practice 2021 0223

import os
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_agg import RendererAgg
_lock = RendererAgg.lock


from skimage import io



from .network_analysis_visualization import (plot_score_per_cluster,
    plot_scores_as_rank, plot_carto, plot_carto_terms, plot_score_comparison)

def network_analysis_01(path_links, embedding_key,
    default_cluster_0, default_cluster_1, cluster_column_name,
    default_gene,
    network_score_kinds, n_genes,
    description_0, title, path_adata=None, lock_matplotlib=True):


    ## Define functions and hyperoarameters

    def first_data_prep_links(path_links):
        if os.path.isdir("tmp"):
            pass
        else:
            os.makedirs("tmp")

        path_score = os.path.join("tmp", os.path.basename(path_links).replace(".celloracle.links", "merged_score.parquet"))
        path_palette = os.path.join("tmp", os.path.basename(path_links).replace(".celloracle.links", "palette.parquet"))

        if os.path.isfile(path_score) & os.path.isfile(path_palette):
            #pass
            print(1)
        else:
            import celloracle as co
            links = co.load_hdf5(path_links)
            links.merged_score.to_parquet(path_score)
            links.palette.to_parquet(path_palette)

    def first_data_prep_adata(path_adata):

        path_emgedding_img = os.path.join("tmp", os.path.basename(path_adata).replace(".h5ad", ".png"))

        if os.path.isfile(path_emgedding_img):
            pass
        else:
            import scanpy as sc
            adata = sc.read_h5ad(path_adata)
            plt.rcParams["savefig.dpi"] = 300
            fig = sc.pl.embedding(adata=adata, basis=embedding_key, color=cluster_column_name,
                                  s=50, return_fig=True, legend_loc="on data")
            fig.savefig(path_emgedding_img)



    @st.cache_resource
    def load_anndata(path_adata):
        #path = "test3.h5ad"
        adata = sc.read_h5ad(path_adata, backed="r")

        return adata

    @st.cache_data
    def load_network_data(path_links):

        path_score = os.path.join("tmp", os.path.basename(path_links).replace(".celloracle.links", "merged_score.parquet"))
        path_palette = os.path.join("tmp", os.path.basename(path_links).replace(".celloracle.links", "palette.parquet"))

        merged_score = pd.read_parquet(path_score)
        palette = pd.read_parquet(path_palette)

        return merged_score, palette

    @st.cache_data
    def load_cluster_image(path_adata):
        path_emgedding_img = os.path.join("tmp", os.path.basename(path_adata).replace(".h5ad", ".png"))
        emgedding_img = io.imread(path_emgedding_img)
        return emgedding_img

    def plot_embeddings(x, args={}):
        fig = sc.pl.embedding(adata=adata, basis=embedding_key, color=x, s=50, return_fig=True, **args)
        return fig




    ### APP starts here

    ## Load data
    first_data_prep_links(path_links)
    merged_score, palette = load_network_data(path_links)
    genes_in_links = merged_score.index.unique()

    first_data_prep_adata(path_adata)
    emgedding_img = load_cluster_image(path_adata)

    ## Side bar

    st.sidebar.write(f"# {title}")

    st.sidebar.write("Please select data below.")

    # Select cluster
    st.sidebar.write("## Cluster")
    clusters = list(palette.index.values)
    cluster = st.sidebar.selectbox("Select cluster to show detailed data",
                                   clusters, index=clusters.index(default_cluster_0))


    st.sidebar.write("## Gene")
    # Select gene
    genes = sorted(np.unique(merged_score.index))
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

        #col1, col2 = st.columns(2)
        st.write("### Clustering")
        st.image(emgedding_img)

        #col2.write("### Gene expression")
        #col2.pyplot(plot_embeddings(gene, args={"use_raw":False, "cmap": "viridis"}))


    st.write(f"# Network scores")
    st.write(f"## 1. Top {n_genes} genes in {cluster} GRN centrality scores")

    col1, col2, col3 = st.columns(3)
    if lock_matplotlib:
        with _lock:
            col1.pyplot(plot_scores_as_rank(merged_score=merged_score,
                                            cluster=cluster,
                                            value="degree_centrality_all",
                                            n_gene=n_genes))
            col2.pyplot(plot_scores_as_rank(merged_score=merged_score,
                                            cluster=cluster,
                                            value="betweenness_centrality",
                                            n_gene=n_genes))
            col3.pyplot(plot_scores_as_rank(merged_score=merged_score,
                                            cluster=cluster,
                                            value="eigenvector_centrality",
                                            n_gene=n_genes))
    else:
        col1.pyplot(plot_scores_as_rank(merged_score=merged_score,
                                            cluster=cluster,
                                            value="degree_centrality_all",
                                            n_gene=n_genes))
        col2.pyplot(plot_scores_as_rank(merged_score=merged_score,
                                        cluster=cluster,
                                        value="betweenness_centrality",
                                        n_gene=n_genes))
        col3.pyplot(plot_scores_as_rank(merged_score=merged_score,
                                        cluster=cluster,
                                        value="eigenvector_centrality",
                                        n_gene=n_genes))


    st.write(f"## 2. All scores in {cluster} GRN")

    df = merged_score[merged_score.cluster == cluster]
    st.dataframe(df)

    st.write(f"## 3. Network score comparison between GRNs")
    st.write(f"### GRN 1: {cluster1},     GRN 2: {cluster2}")
    st.write(f"### Network score : {value}")
    st.plotly_chart(plot_score_comparison(merged_score=merged_score,
                                          value=value,
                                          cluster1=cluster1,
                                          cluster2=cluster2))
    
    st.write("# Network dynamics")
    st.write(f"### Selected gene: {gene}")
    if gene in genes_in_links:
        pass
    else:
        st.write(f"### Caution! {gene} does not have any connection in the filtered GRNs.")
    st.write(f"## 1. Score dynamics")
    if lock_matplotlib:
        with _lock:
            st.pyplot(plot_score_per_cluster(merged_score=merged_score, 
                                             palette=palette, gene=gene,
                                             figsize=[6, 5]))
    else:
        st.pyplot(plot_score_per_cluster(merged_score=merged_score, 
                                             palette=palette, gene=gene,
                                             figsize=[6, 5]))
        
    st.write("## 2. Gene cartography analysis")
    col1, col2 = st.columns([1, 2])
    col1.write("### 2.1 Cartography summary")
    if gene in genes_in_links:
        if lock_matplotlib:
            with _lock:
                col1.pyplot(plot_carto_terms(merged_score=merged_score, 
                                             palette=palette, gene=gene))
        else:
            col1.pyplot(plot_carto_terms(merged_score=merged_score, 
                                             palette=palette, gene=gene))
    else:
        col1.write(f"{gene} does not have any connection in the filtered GRNs.")
    col2.write(f"### 2.2 Cartography of {gene} in {cluster} GRN")
    if lock_matplotlib:
        with _lock:
            col2.pyplot(plot_carto(merged_score=merged_score,
                                   gene=gene, cluster=cluster))
    else:
        col2.pyplot(plot_carto(merged_score=merged_score,
                                   gene=gene, cluster=cluster))


