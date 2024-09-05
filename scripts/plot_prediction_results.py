import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import AutoMinorLocator


def merge_exp_pred_matrix_models(pred_prot, matrix_type):
    """
    Merge DataFrame with experimental abundance scores and
    DataFrame with abundance scores predicted using each
    of the four substitution matrix models. 

    pred_prot : name of protein for which to perform predictions
    matrix_type : specify whether to use mean or median abundance score.
    """
    # read vampseq data and get data for pred_prot specifically
    vamp_df = pd.read_csv("../data/vampseq_data/processed/vampseq_soluble_df.csv")
    pred_prot_vamp_df = vamp_df.groupby("prot").get_group(pred_prot)

    # load predictions 
    pred_df = pd.read_csv("../output/predictions_from_matrices/dimers/"+
                          f"pred_{pred_prot}_from_leave_out_{pred_prot}_{matrix_type}.csv")

    # merge data 
    pred_df = pred_prot_vamp_df[["variant","abundance_score"]].merge(pred_df, how="left", left_on="variant", right_on="variant")
    
    return pred_df


def make_eval_dfs(matrix_type):
    """
    Calculate Pearson correlation and MAE between experimental and predicted abundance 
    scores for all six proteins and all four substitution matrix-based models and collect
    all model evaluations in new DataFrames that can be used for plotting. 
    """
    # make empty DataFrames to store results
    df_corr = pd.DataFrame()
    df_MAE = pd.DataFrame()

    for pred_prot in ['PTEN','TPMT','CYP2C9','NUDT15','ASPA','PRKN']:

        # get DataFrame with experimental and predicted abundance scores
        pred_df = merge_exp_pred_matrix_models(pred_prot, matrix_type)

        # calculate pearson correlation between experimental and 
        # predicted abundance scores for all four substitution matrix models
        corr_series = pred_df.corr(method="pearson")["abundance_score"][1:]

        # calculate MAE between experimental and 
        # predicted abundance scores for all four substitution matrix models
        MAE_series = np.abs(pred_df.drop(columns=["variant"]).sub(pred_df["abundance_score"],axis=0)).mean(axis=0)[1:]

        # add results to DataFrames
        df_corr[pred_prot] = corr_series
        df_MAE[pred_prot] = MAE_series

    return df_corr, df_MAE


def corr_df_training_combinations(pred_prot,struc_type,matrix_type): 
    """
    Make DataFrame that for each combination of training data proteins contains 
    the correlation between pred_prot abundance data and the pred_prot predictions 
    obtained from each training data combination.  
    struc_type : either exposure or combined
    matrix_type : either mean or median
    """
    # read vampseq data and get data for pred_prot specifically
    vamp_df = pd.read_csv("../data/vampseq_data/processed/vampseq_soluble_df.csv")
    pred_prot_vamp_df = vamp_df.groupby("prot").get_group(pred_prot)

    # read predictions for pred_prot made with increasing amounts of data
    combinations_pred_df = pd.read_csv("../output/predictions_from_matrices/dimers_combinations/"+
                                       f"pred_{pred_prot}_{struc_type}_combinations_{matrix_type}.csv",index_col=0)

    # merge experimental data and predictions
    merged_df = pred_prot_vamp_df[["prot","variant","abundance_score"]].merge(combinations_pred_df,how="left",
                                                                              left_on="variant",right_on="variant")

    # calculate correlation between abundance and predictions from all training data combinations
    corr_series = merged_df.corr(method="pearson")['abundance_score']

    # make dataframe with correlations for plotting
    number_of_datasets = [int(idx[0]) for idx in corr_series.index[1:]]
    corr_df = pd.DataFrame(corr_series[1:])
    corr_df["number_of_datasets"] = number_of_datasets 
    corr_df = corr_df.rename(columns={"abundance_score":"Pearson correlation"})

    # calculate mean correlation for every number of datasets 
    corr_mean_df = corr_df.groupby("number_of_datasets").mean().reset_index()

    return corr_df, corr_mean_df
    
    
def plot_eval_dfs(eval_df, eval_metric):
    
    eval_df = eval_df.rename(index={'global':"Global",'exposure':"Exposure",
                                    'secondary_structure':"Secondary",'combined':"Combined"})
    
    eval_df = eval_df.stack().reset_index()
    eval_df.columns = ["Model","Protein",eval_metric]

    # set seaborn theme
    sns.set_theme(style="ticks", palette="colorblind", font_scale=1.5)

    # make swarmplot with evaluation metrics for each model and protein
    plt.figure(figsize=(5.5,8.5))
    swarm = sns.swarmplot(data=eval_df,x="Model",y=eval_metric,hue="Protein",
                      size=12,linewidth=2,edgecolor="gray",palette="colorblind")

    # also show medians for each model in plot
    sns.boxplot(data=eval_df,x="Model",y=eval_metric,
                showmeans=False,meanline=False,
                medianprops={'visible': True, 'color': 'k', 'ls': '-', 'lw': 3},
                whiskerprops={'visible': False},
                showfliers=False,showbox=False,
                showcaps=False,
                ax=swarm)

    # define position of legend
    #plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0,edgecolor='black')
    plt.legend(ncol=2,loc='upper center', bbox_to_anchor=(0.5, 1.3))
    
    # draw vertical lines to seperate results from different models
    minor_locator = AutoMinorLocator(2)
    plt.gca().xaxis.set_minor_locator(minor_locator)
    plt.grid(which='minor')
    
    plt.xticks(rotation=-45)
    plt.ylim(0.2,0.65)
    
    plt.tight_layout()
    plt.savefig(f"../figures/predictions_from_matrices/dimers/eval_all_{eval_metric}_{matrix_type}.pdf")
    plt.close()

    
def plot_corr_hexbin(struc_type, matrix_type):
    """
    Make scatter plots with experimental vs. predicted abundance scores 
    for all six proteins. The type of matrices used for predictions are 
    defined by struc_type.
    struc_type : global, secondary_structure, exposure or combined.
    matrix_type : mean or median
    """
    proteins = ['PTEN','TPMT','CYP2C9','NUDT15','ASPA','PRKN']
    
    plt.rcParams.update({'font.size': 13})
    plt.figure(figsize=(11,7))
    fs = 13
    i = 0

    for pred_prot in proteins:

        # get DataFrame with experimental and predicted abundance scores
        pred_df = merge_exp_pred_matrix_models(pred_prot, matrix_type)
        
        # filter to keep only predictions for struc_type
        pred_df = pred_df[["abundance_score",struc_type]].dropna()
        
        # calculate pearson between exp and pred
        corr_r = np.around(pred_df.corr(method="pearson")["abundance_score"][struc_type],2)
        corr_rs = np.around(pred_df.corr(method="spearman")["abundance_score"][struc_type],2)
        
        x = pred_df["abundance_score"].values
        y = pred_df[struc_type].values

        i += 1
        ax = plt.subplot(2, 3, i)  
        plt.hexbin(x,y,mincnt=1,gridsize=60,bins='log',
                   vmin=1,vmax=50,cmap='Spectral_r',
                   label=r'$\rho_r$ = '+str(corr_r))
        plt.title(pred_prot,fontsize=14)
        plt.xlabel('Experimental abundance',fontsize=fs)
        plt.ylabel('Predicted abundance',fontsize=fs)
        plt.ylim(plt.xlim())
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5)) 
        #plt.colorbar()
        #plt.legend(frameon=False,handlelength=0,handleheight=0)
        plt.annotate(f"$r$ = {str(corr_r)}, $r_s$ = {str(corr_rs)}", xy=(0.075,0.87), xycoords='axes fraction', fontsize=14)
        
    plt.tight_layout()
    plt.savefig(f"../figures/predictions_from_matrices/dimers/exp_pred_scatter_{struc_type}_{matrix_type}.pdf")
    plt.close()    
    

def plot_corr_combinations(struc_type, matrix_type):

    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(12,7))
    c = 0

    for pred_prot in ['PTEN','TPMT','CYP2C9','NUDT15','ASPA','PRKN']:
        c += 1
        plt.subplot(2,3,c)

        corr_df, corr_mean_df = corr_df_training_combinations(pred_prot,struc_type,matrix_type)

        plt.plot(corr_mean_df['number_of_datasets'] - 1, corr_mean_df['Pearson correlation'],
                 lw=0.5, marker='o', color='black', ms=13)

        sns.swarmplot(data=corr_df, x='number_of_datasets', y='Pearson correlation',
                      s=6.5, linewidth=1, edgecolor='gray',
                      color=sns.color_palette()[1])
        
        if struc_type == "exposure":
            if pred_prot == 'ASPA':
                plt.ylim(0.3,0.49)
            elif pred_prot == 'PRKN':
                plt.ylim(0.37,0.56)
            else:
                plt.ylim(0.4,0.59)

        if struc_type == "combined":
            if pred_prot == "PRKN":
                plt.ylim(0.25,0.55)
            else:
                plt.ylim(0.3,0.6)
                
        plt.xlabel('Number of datasets')
        plt.ylabel(r'$r$')
        plt.title(pred_prot,fontsize=14)

    plt.tight_layout()
    plt.savefig("../figures/predictions_from_matrices/dimers_combinations/"+
                f"corr_combinations_{struc_type}_{matrix_type}.pdf")    
     
    
def make_grid_result_matrix(grid_vamp_series, wcn_cut_list, rasa_cut_list):
    """
    Convert Pandas Series grid_vamp_series with correlation coefficients or 
    MAEs calculated for all wcn_cut and rasa_cut combinations to Numpy matrix. 
    """
    matrix = np.zeros((len(wcn_cut_list),len(rasa_cut_list)))
    
    for i in range(len(wcn_cut_list)):
        
        for j in range(len(rasa_cut_list)):
            
            wcn_cut = wcn_cut_list[i]
            rasa_cut = rasa_cut_list[j]
            
            matrix[i][j] = grid_vamp_series[f"{wcn_cut}_{rasa_cut}"]
            
    return matrix


def analyse_gridsearch_predictions(pred_prot, matrix_type, wcn_cut_list, rasa_cut_list):
    """
    Calculate pearson's correlation coefficient and MAE between experimental and 
    predicted abundance scores for every wcn_cut and rasa_cut combination. 
    """
    # load gridsearch predictions for pred_prot
    grid_df = pd.read_csv("../output/predictions_from_matrices/dimers_gridsearch/"+
                          f"pred_{pred_prot}_from_leave_out_{pred_prot}_{matrix_type}.csv")

    # load vampseq data and group by prot
    vamp_df = pd.read_csv("../data/vampseq_data/processed/vampseq_soluble_df.csv")
    pred_prot_vamp_df = vamp_df.groupby("prot").get_group(pred_prot)

    # merge gridsearch predictions and vampseq data 
    grid_vamp_df = pred_prot_vamp_df[["variant","abundance_score"]].merge(grid_df,how="left",on=["variant"])

    # calculate pearson's correlation coefficient between experimental and predicted scores for all cutoff combinations
    grid_vamp_corr = grid_vamp_df.corr(method='pearson')['abundance_score']

    # calculate MAE between experimental and predicted scores for all cutoff combinations
    grid_vamp_mae = np.abs(grid_vamp_df.drop(columns=["variant"]).sub(grid_vamp_df["abundance_score"],axis=0)).mean(axis=0)
    
    # put results into matrix format
    corr_matrix = make_grid_result_matrix(grid_vamp_corr, wcn_cut_list, rasa_cut_list)
    mae_matrix = make_grid_result_matrix(grid_vamp_mae, wcn_cut_list, rasa_cut_list)
    
    return corr_matrix, mae_matrix
    
    
def plot_gridsearch_heatmap(matrix, rasa_cut_list, wcn_cut_list, plot_title, figname, vmin, vmax, cmap):
    
    plt.rcParams['font.size'] = '16'
    
    fig, ax = plt.subplots(figsize=(8,5))
    im = ax.imshow(matrix,origin='lower',vmin=vmin,vmax=vmax,cmap=cmap)

    ax.set_xticks(range(len(rasa_cut_list)), labels=rasa_cut_list)
    ax.set_yticks(range(len(wcn_cut_list)), labels=wcn_cut_list)

    ax.set_xlabel("rASA cutoff")
    ax.set_ylabel("WCN cutoff")
    
    ax.set_title(plot_title)
    
    fig.colorbar(im, orientation='vertical',label=plot_title)
    
    for i in range(len(wcn_cut_list)):
        for j in range(len(rasa_cut_list)):
            text = ax.text(j, i, np.around(matrix[i, j],3),fontsize=12,
                           ha="center", va="center", color="w") 
    
    plt.tight_layout()
    plt.savefig(f"../figures/predictions_from_matrices/dimers_gridsearch/{figname}.pdf")
    plt.close()    
    
    
def plot_all_gridseach_results(matrix_type):    

    wcn_cut_list = [0,5,10,12.5,15,17.5,20]
    rasa_cut_list = [0,0.1,0.2,0.3,0.4,0.5]

    corr_matrix_list = []
    mae_matrix_list = []

    for pred_prot in ["PTEN","TPMT","CYP2C9","NUDT15","ASPA","PRKN"]:

        # get results in matrix format
        corr_matrix, mae_matrix = analyse_gridsearch_predictions(pred_prot, matrix_type, wcn_cut_list, rasa_cut_list)

        # plot results per protein
        #plot_gridsearch_heatmap(corr_matrix, rasa_cut_list, wcn_cut_list, 
        #                        f"corr_matrix_{pred_prot}_{matrix_type}", 
        #                        f"corr_matrix_{pred_prot}_{matrix_type}",'viridis')
        #plot_gridsearch_heatmap(mae_matrix, rasa_cut_list, wcn_cut_list, 
        #                        f"MAE_matrix_{pred_prot}_{matrix_type}", 
        #                        f"MAE_matrix_{pred_prot}_{matrix_type}",'viridis')

        # append to lists
        corr_matrix_list.append(corr_matrix)
        mae_matrix_list.append(mae_matrix)

    # get mean correlation and MAE across all proteins    
    corr_matrix_mean = np.mean(corr_matrix_list,axis=0)
    mae_matrix_mean = np.mean(mae_matrix_list,axis=0)

    # plot mean results 
    plot_gridsearch_heatmap(corr_matrix_mean, rasa_cut_list, wcn_cut_list, 
                            fr"$r$", f"corr_matrix_all_{matrix_type}",
                            0.33,0.52,'viridis_r')
    plot_gridsearch_heatmap(mae_matrix_mean, rasa_cut_list, wcn_cut_list, 
                            f"Mean absolute error", f"MAE_matrix_all_{matrix_type}",
                            0.26,0.295,'viridis')    

    
if __name__ == "__main__":    
    
    for matrix_type in ["mean","median"]:
        
        # make swarm plots with correlation coefficients and MAEs for all models and proteins 
        df_corr, df_MAE = make_eval_dfs(matrix_type)
        plot_eval_dfs(df_MAE, "MAE")
        plot_eval_dfs(df_corr, r"$r$")
        
        # make scatter plot with exp vs pred abundance from substitution matrices
        for struc_type in ["global","secondary_structure","exposure","combined"]:
            plot_corr_hexbin(struc_type, matrix_type)
            
        # set seaborn theme
        sns.set_theme(style="ticks", palette="Set2", font_scale=1.3)    
            
        # plot correlation as function of number of datasets
        for struc_type in ["exposure","combined"]:      
            plot_corr_combinations(struc_type, matrix_type)
            
        # plot gridsearch results 
        plot_all_gridseach_results(matrix_type)