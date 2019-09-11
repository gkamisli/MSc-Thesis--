def plot_attention_matrix(attention_weights, figsize, xticklabels):
    
    """
    Function to plot heatmap of attention matrices 
    
    Arguments:
    attention_weights -- a Pandas DataFrame for attention weight matrix, shape of (time_steps, num_examples)
    figsize -- figure size defined in the seaborn.set()
    xticklabels -- Show xticklabels, if True; otherwise, do not show xticklabels
    
    Returns:
    Heatmap plot
    """
    
    import seaborn as sns
    import numpy as np
    import pandas as pd
    sns.set(rc={'figure.figsize':(15,15)})
    
    sns.set(rc={'figure.figsize':figsize})
    sns.heatmap(attention_weights, xticklabels=xticklabels, yticklabels=True, 
                cbar_kws = dict(use_gridspec=False,location="right", shrink=0.5))