import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import itertools


def plot_substitution_matrix_heatmap(sub_matrix, vmin, vmax, figname, cmap='rocket'):
    """
    Plot substitution matrix as heatmap. 
    sub_matrix : substitution matrix in DataFrame format. 
    """
    plt.figure(figsize=(8,8))
    sns.set(font_scale=1.5)
    m = sns.heatmap(sub_matrix,square=True,vmin=vmin, vmax=vmax, cmap=cmap,
                    cbar_kws={'label': 'Average abundance score','shrink': 0.73}) 
    m.set_facecolor('silver')
    plt.xlabel('Mutation to')
    plt.ylabel('Mutation from')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig('../figures/substitution_matrices/'+figname+'.pdf') 
    plt.close()
    
    
def plot_substitution_matrix_counts(sub_matrix_count,figname):
    """
    Make barplot of the average number of variants contributing 
    to each entry in a "mutation from" row in a substitution matrix. 
    sub_matrix_count : DataFrame with variant counts for substitution matrix.
    """
    # set synonymous scores to nan
    for aa1,aa2 in zip(sub_matrix_count.index,sub_matrix_count.index):
        sub_matrix_count[aa1][aa2] = np.nan

    # calculate average counts for matrix rows
    mean_count_series = sub_matrix_count.mean(axis=1,skipna=True)

    # make barplot 
    sns.reset_orig() 
    fig = plt.figure(figsize=(3,8))
    plt.rcParams.update({'font.size': 19})
    plt.barh(np.flip(mean_count_series.index),np.flip(mean_count_series),color='black')
    plt.xlabel('Average number \n of variants')
    plt.locator_params(axis='x', nbins=5)
    plt.margins(y=0)
    plt.tight_layout()
    plt.savefig('../figures/substitution_matrices/'+figname+'.pdf') 
    plt.close()
    
    
def plot_substitution_matrix_symmetry(df,vmin,vmax,filename):
    """
    Plot heatmap showing abundance score average differences
    between symmetric substitutions. 
    df : DataFrame with average mutational scores.
    """
    # calculate symmetry data
    diff_dict = {}
    for amino_acid in df.index:
        diff_dict[amino_acid] = {}

    for aa1,aa2 in itertools.product(df.index,repeat=2):
        diff = df[aa1][aa2] - df[aa2][aa1] # mutation to minus mutation from
        diff_dict[aa1][aa2] = diff

    symmetry_df = pd.DataFrame(diff_dict).reindex(index=df.index,columns=df.columns)

    # plot symmetry data
    plt.figure(figsize=(8,8))
    sns.set(font_scale=1.5)
    m = sns.heatmap(symmetry_df,square=True,vmin=vmin,vmax=vmax,cmap='bwr',
                cbar_kws={'label': 'Average abundance score difference','shrink': 0.73}) 
    m.set_facecolor('dimgray')
    plt.xlabel('Mutation to')
    plt.ylabel('Mutation from')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig('../figures/substitution_matrices/'+filename+'.pdf') 
    plt.close()
    

def plot_matrices_dimers(matrix_name, matrix_type):
    """
    Plot subsitution matrices heatmaps, barplots with variant counts,
    and substitution matrix symmetry heatmaps for all structure types
    (Global, Exposure, Secondary, Combined). Plots are made for
    substitution matrices constructed with data from all six VAMP-seq
    datasets and using ASPA and NUDT15 dimer structures for structure
    feature calculation.

    matrix_name : "all" 
    matrix_type : mean or median
    """
    # get global matrix with averages across all proteins    
    global_mean_matrix = pd.read_pickle("../output/substitution_matrices/dimers/global/all/mean.pkl")    

    # calculate mean score per row to define aa order in all other plots
    mut_from_mean = np.nanmean(global_mean_matrix,axis=1)
    global_aa_sort = np.flip(global_mean_matrix.index[np.argsort(mut_from_mean)])       

    # define structural environments for which to plot data
    struc_types_list = ["buried","buried_helix","buried_loop","buried_strand","exposed", "exposed_helix",
                        "exposed_loop","exposed_strand","global","helix","loop","strand"]

    # make plots for all struc_types
    for struc_type in struc_types_list:

        # create output directory 
        if not os.path.isdir(f"../figures/substitution_matrices/dimers/{struc_type}"):
            os.mkdir(f"../figures/substitution_matrices/dimers/{struc_type}")       
        if not os.path.isdir(f"../figures/substitution_matrices/dimers/{struc_type}/{matrix_name}"):
            os.mkdir(f"../figures/substitution_matrices/dimers/{struc_type}/{matrix_name}")       
        
        # get substitution matrix with average abundance scores
        matrix = pd.read_pickle(f"../output/substitution_matrices/dimers/{struc_type}/{matrix_name}/{matrix_type}.pkl")

        # sort residues according to global_aa_sort scale before plotting 
        matrix_sorted = matrix[np.flip(global_aa_sort)].loc[global_aa_sort]

        # plot sorted matrix as heatmap
        plot_substitution_matrix_heatmap(matrix_sorted, 0, 1, f"dimers/{struc_type}/{matrix_name}/{matrix_type}", cmap='rocket')

        # plot symmetry heatmap
        plot_substitution_matrix_symmetry(matrix_sorted,-0.5,0.5,f"dimers/{struc_type}/{matrix_name}/{matrix_type}_symmetry")
        
        # get matrix with subsitution type counts
        matrix_counts = pd.read_pickle(f"../output/substitution_matrices/dimers/{struc_type}/{matrix_name}/count.pkl") 

        # sort residues according to global_aa_sort scale before plotting 
        matrix_counts_sorted = matrix_counts[np.flip(global_aa_sort)].loc[global_aa_sort]

        # make barplot with average variant counts per row in matrix
        plot_substitution_matrix_counts(matrix_counts_sorted,f"dimers/{struc_type}/{matrix_name}/{matrix_type}_count")
        
        
if __name__ == "__main__":

    for matrix_type in ["mean","median"]:
        
        # make plots for matrices calculated with all data
        # and for matrices calculated without ASPA and without PRKN data
        for matrix_name in ["all","leave_out_ASPA","leave_out_PRKN"]: 
        
            plot_matrices_dimers(matrix_name, matrix_type)
