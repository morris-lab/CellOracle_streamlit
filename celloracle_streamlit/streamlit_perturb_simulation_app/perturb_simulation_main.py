# Streamlit practice 2021 0223

import streamlit as st
import pandas as pd
import numpy as np
import scanpy as sc
import os


from ..applications import Oracle_development_module, Oracle_systematic_analysis_helper


from .perturb_simulation_visualization import (plot_cluster_and_dev_flow,
                                           plot_sim_vectorfield,
                                           plot_sim_quiver,
                                           plot_ip_score,
                                           plot_ip_distribution,
                                           plot_selected_pseudotime)

def perturb_simulation_set_01(path_adata, path_sim_data, embedding_key, cluster_column_name,
    default_gene, default_unit,
    dev_scale_min, dev_scale_max, dev_scale_step, dev_scale_default,
    sim_scale_min, sim_scale_max, sim_scale_step, sim_scale_default, sim_coef,
    ip_min, ip_min_default, ip_max, ip_max_default, ip_step,
    description_0, description_1, title,
    pseudotime_min=None, pseudotime_min_default=None, pseudotime_max=None, pseudotime_max_default=None):


    def first_data_prep_ip(path_sim_data):
        if os.path.isdir("tmp"):
            pass
        else:
            os.makedirs("tmp")

        path_nip = os.path.join("tmp", os.path.basename(path_sim_data).replace(".hdf5", "negative_ip_sum.parquet"))

        if os.path.isfile(path_nip):
            pass
        else:
            helper = Oracle_systematic_analysis_helper(hdf5_file_path=path_sim_data)
            helper.get_negative_ip_sum_for_all_data(verbose=False, return_result=False)
            helper.negative_ip_sum.to_parquet(path_nip)


    ## Define functions and hyperoarameters
    @st.cache(allow_output_mutation=True)
    def load_data(path_adata, path_sim_data):
        #path = "test3.h5ad"
        adata = sc.read_h5ad(path_adata)

        if adata.raw is not None:
            adata.raw = None
            adata.write_h5ad(path_adata)
            st.write("adata raw overwriting")
        if len(adata.layers.keys()) != 0:
            adata.layers = None
            adata.write_h5ad(path_adata)
            st.write("adata layers overwriting")

        # Make a new Oracle_developmennt_module object to read calculation result
        dev = Oracle_development_module()
        dev.set_hdf_path(path=path_sim_data)
        meta_data = dev.get_hdf5_info()

        return adata, dev, meta_data


    @st.cache(allow_output_mutation=True, suppress_st_warning=True)
    def load_ip_scores(path_sim_data):
        # Load data with Oracle_systematic_analysis_helper.
        path_nip = os.path.join("tmp", os.path.basename(path_sim_data).replace(".hdf5", "negative_ip_sum.parquet"))
        helper = Oracle_systematic_analysis_helper(hdf5_file_path=path_sim_data)
        helper.negative_ip_sum = pd.read_parquet(path_nip)
        return helper

    #@st.cache(allow_output_mutation=True)
    def load_oravle_dev(gene, misc):
        dev.load_hdf5(gene=gene, misc=misc)
        return dev

    def plot_embeddings(x, args={}):
        fig = sc.pl.embedding(adata=adata, basis=embedding_key,
                              color=x, s=50, return_fig=True, **args)
        return fig


    ### APP starts here

    ## Load data
    adata, dev, meta_data= load_data(path_adata=path_adata, path_sim_data=path_sim_data)

    if pseudotime_min is not None:
        first_data_prep_ip(path_sim_data=path_sim_data)
        helper = load_ip_scores(path_sim_data=path_sim_data)

    ## Side bar

    st.sidebar.write(f"# {title}")

    st.sidebar.write("## 1. Check gene expression")
    # Select gene
    genes = sorted(adata.var.index.values)
    gene_viz = st.sidebar.selectbox("Select gene to show gene expression",
                                genes,
                                index=genes.index(default_gene))

    st.sidebar.write("## 2. Simulation")

    # Select metric
    units = list(meta_data["misc_list"])
    unit = st.sidebar.selectbox("Select simulation unit",
                                 units,
                                 index=units.index(default_unit))

    #st.sidebar.write("Select gene from the list")
    genes_sim = sorted(meta_data["misc_gene_dictionary"][unit])
    gene = st.sidebar.selectbox("Select gene for simulation",
                                genes_sim,
                                index=genes_sim.index(default_gene))

    dev = load_oravle_dev(gene, unit)

    # Parameters for visualization
    scale_devflow = st.sidebar.slider("Scale parameter for development flow vectors",
                                      min_value=dev_scale_min, max_value=dev_scale_max,
                                      step=dev_scale_step, value=dev_scale_default)

    scale_perturb = st.sidebar.slider("Scale parameter for perturb simulation vectors",
                                      min_value=sim_scale_min, max_value=sim_scale_max,
                                      step=sim_scale_step, value=sim_scale_default)

    vmin = st.sidebar.slider("Perturbation score min",
                           min_value=ip_min,
                           max_value=-ip_step, step=ip_step, value=ip_min_default)
    vmax = st.sidebar.slider("Perturbation score max",
                           min_value=ip_step,
                           max_value=ip_max, step=ip_step, value=ip_max_default)


    if pseudotime_min is not None:

        st.sidebar.write("## 3. Sort genes by Negative PS sum")
        st.sidebar.write("Select digitized pseudotime range for PS sum calculation")
        lower = st.sidebar.number_input("Lower",
                                        min_value=pseudotime_min,
                                        max_value=pseudotime_max,
                                        step=1,
                                        value=pseudotime_min_default)

        upper = st.sidebar.number_input("Upper",
                                        min_value=pseudotime_min,
                                        max_value=pseudotime_max,
                                        step=1,
                                        value=pseudotime_max_default)

    ## Main
    st.write(f"# {title}")
    st.write("## About")
    st.write(description_0)

    st.write("# 1. scRNA-seq data expolation")

    col1, col2 = st.columns(2)
    col1.write("### Cell type annotation")
    col1.pyplot(plot_embeddings(cluster_column_name, args={"legend_loc": "on data"}))

    col2.write("### Gene expression")
    col2.pyplot(plot_embeddings(gene_viz, args={"use_raw":False, "cmap": "viridis"}))

    st.write(f"# 2. CellOracle gene perturbation simulation results")
    st.write(description_1)

    st.write(f"### Selected Lineage for simulation: {unit}, Selected Gene: {gene}")

    st.write(f"## 1. Development direction")
    st.pyplot(plot_cluster_and_dev_flow(dev, scale_devflow))

    st.write(f"## 2. KO Simulation results")
    st.pyplot(plot_sim_quiver(dev, gene, scale_perturb, coef=sim_coef))
    st.pyplot(plot_sim_vectorfield(dev, gene, scale_perturb))

    st.write("## 3. Perturbation score (PS)")
    st.pyplot(plot_ip_score(dev, vmin, vmax, scale_perturb))
    st.pyplot(plot_ip_distribution(dev, vmin, vmax))

    if pseudotime_min is not None:

        st.write(f"# 3. Sort genes by Negative PS sum")

        range_int = list(np.arange(int(lower), int(upper) + 1))
        range_ = [str(i) for i in range_int]
        range_ = ",".join(range_)

        st.write("Selected simulation unit: ", unit)
        st.write("Selected pseudotime range: " + range_)
        col1, col2 =  st.columns([1.2, 1])
        col1.write("Gene list sorted by sum of PS")
        col1.dataframe(helper.sort_TFs_by_neagative_ip(unit, pseudotime=range_))
        col2.write("Selected grid points")
        col2.pyplot(plot_selected_pseudotime(dev, pseudotime_selected=range_int))
