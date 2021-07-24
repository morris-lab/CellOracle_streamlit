# Custom wrapper function for streamlit visualization
import matplotlib.pyplot as plt

FIG_SIZE_1 = [5, 5]
FIG_SIZE_2 = [10, 5]
S=5
S_GRID=50


def plot_cluster_and_dev_flow(dev, scale_devflow):


    fig, ax = plt.subplots(1, 2, figsize=FIG_SIZE_2)

    dev.plot_cluster_cells_use(ax=ax[0], s=S)
    ax[0].set_title("Clustering")
    dev.plot_reference_flow_on_grid(ax=ax[1], scale=scale_devflow, s=S)
    ax[1].set_title("Vectorfield: Development pseudotime")

    return fig

def plot_sim_vectorfield(dev, gene, scale_perturb):

    fig, ax = plt.subplots(1, 2, figsize=FIG_SIZE_2)

    dev.plot_simulation_flow_on_grid(ax=ax[0], scale=scale_perturb, s=S)
    ax[0].set_title(f"Vectorfield: {gene} KO simulation")
    dev.plot_simulation_flow_random_on_grid(ax=ax[1], scale=scale_perturb, s=S)
    ax[1].set_title(f"Vectorfield: Simulation with randomized GRNs")

    return fig

def plot_sim_quiver(dev, gene, scale_perturb, coef=40):

    fig, ax = plt.subplots(1, 2, figsize=FIG_SIZE_2)

    dev.plot_quiver(ax=ax[0], scale=scale_perturb*coef, s=S)
    ax[0].set_title(f"Quiver plot: {gene} KO simulation")
    dev.plot_quiver_random(ax=ax[1], scale=scale_perturb*coef, s=S)
    ax[1].set_title(f"Quiver plot: Simulation with randomized GRNs")

    return fig

def plot_ip_score(dev, vmin, vmax, scale_perturb):

    fig, ax = plt.subplots(1, 2, figsize=FIG_SIZE_2)
    dev.plot_inner_product_on_grid(ax=ax[0], vmin=vmin, vmax=vmax)
    ax[0].set_title("Inner-product score")

    dev.plot_inner_product_on_grid(ax=ax[1], vmin=vmin, vmax=vmax)
    dev.plot_simulation_flow_on_grid(ax=ax[1], scale=scale_perturb, show_background=False)
    ax[1].set_title("Inner-product score \nwith KO simulation vectorfield")

    return fig

def plot_ip_distribution(dev, vmin, vmax):
    fig, ax = plt.subplots(1, 2, figsize=FIG_SIZE_2)

    dev.plot_inner_product_on_pseudotime(ax=ax[0], vmin=vmin, vmax=vmax, s=S_GRID)

    dev.plot_inner_product_as_box(ax=ax[1], vmin=vmin, vmax=vmax)

    return fig

def plot_selected_pseudotime(dev, pseudotime_selected=[]):
    fig, ax = plt.subplots(figsize=FIG_SIZE_1)
    dev.plot_selected_pseudotime_on_grid(ax=ax, pseudotime_selected=pseudotime_selected)
    return fig
