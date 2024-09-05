import numpy as np
import pandas as pd
import itertools
import pickle
import os


def make_empty_matrix(dim):
    
    matrix = np.zeros((dim,dim))
    matrix[:] = np.nan
    
    return matrix


def calc_substitution_matrix(df,outdir):
    """
    Given input DataFrame with abundance_score, aa_ref and aa_var columns, 
    calculate mean and median abundance score for all 380 substitution types. 
    """
    # check that input data doesn't contain nan entries
    assert ~np.any(np.isnan(df.abundance_score))

    # define amino acids to base calculation on
    amino_acid_list = ['A','C','D','E','F','G','H','I','K','L',
                       'M','N','P','Q','R','S','T','V','W','Y']
    aa_len = len(amino_acid_list)
    
    # make empty matrices
    mean_matrix = make_empty_matrix(aa_len)
    median_matrix = make_empty_matrix(aa_len)
    count_matrix = make_empty_matrix(aa_len)
    
    # make arrays with relevant data
    aa_ref_arr = df.aa_ref.values
    aa_var_arr = df.aa_var.values
    abundance_score_arr = df.abundance_score.values

    # calculate mean, median and number of scores for all aa combinations
    for i,j in itertools.product(range(aa_len),range(aa_len)):

        aa_ref = amino_acid_list[i]
        aa_var = amino_acid_list[j]

        score_arr = abundance_score_arr[(aa_ref_arr == aa_ref) & (aa_var_arr == aa_var)] 

        abundance_count = len(score_arr)

        if abundance_count > 0:
            abundance_mean = np.mean(score_arr)
            abundance_median = np.median(score_arr)

        else:
            abundance_mean = np.nan
            abundance_median = np.nan

        mean_matrix[i][j] = abundance_mean
        median_matrix[i][j] = abundance_median
        count_matrix[i][j] = abundance_count

    # create output directory 
    if not os.path.isdir(f"../output/substitution_matrices/{outdir}"):
        os.mkdir(f"../output/substitution_matrices/{outdir}")       
    
    # save results
    for matrix_name,matrix in zip(["mean","median","count"],[mean_matrix,median_matrix,count_matrix]):
        
        # make df from matrix, aa_ref is index, aa_var is columns
        matrix_to_df = pd.DataFrame(matrix,columns=amino_acid_list,index=amino_acid_list)
        matrix_to_df.to_pickle(f"../output/substitution_matrices/{outdir}/{matrix_name}.pkl")
        
        
def calc_matrices_for_selection(vamp_struc_df, selection_arr, outdir_pre):
    """
    Use calc_substitution_matrix to calculate substitution matrices using 
    abundance scores in vamp_struc_df for all proteins combined, individual proteins, 
    and all leave-out-one protein combinations. Use structural features given
    by the boolean array selection_arr, which specifices for each entry in vamp_struc_df 
    if this variant should be used for the matrix calculation or not. 
    """
    # define proteins
    proteins = ['PTEN','TPMT','CYP2C9','NUDT15','ASPA','PRKN']
    
    # get protein name array from dataframe 
    prot_arr = vamp_struc_df.prot.values 
    
    # calculate matrix using data from all proteins
    vamp_struc_df_filtered = vamp_struc_df[selection_arr] 
    calc_substitution_matrix(vamp_struc_df_filtered, f"{outdir_pre}/all")

    for prot in proteins:

        # calculate matrix using data for single proteins
        prot_bool_arr = np.isin(prot_arr, [prot])
        vamp_struc_df_filter = vamp_struc_df[prot_bool_arr & selection_arr]
        calc_substitution_matrix(vamp_struc_df_filter,f"{outdir_pre}/{prot}")

        # calculate global matrices leaving out single proteins
        proteins_leave_one_out = proteins.copy()
        proteins_leave_one_out.remove(prot)
        prot_bool_arr = np.isin(prot_arr, [proteins_leave_one_out])
        vamp_struc_df_filter = vamp_struc_df[prot_bool_arr & selection_arr]
        calc_substitution_matrix(vamp_struc_df_filter,f"{outdir_pre}/leave_out_{prot}")        
    
        
def calc_matrices_for_selection_all_combinations(vamp_struc_df, selection_arr, combinations, outdir_pre):
    """
    Calculate substitution matrices using different combinations of
    VAMP-seq datasets, as specificed in the combinations dictionary. 
    """
    # get protein name array from dataframe 
    prot_arr = vamp_struc_df.prot.values
    
    # loop through protein combinations, calc matrix for each combination
    for i in combinations:

        # calc matrix using data for all proteins in combinations[i]
        prot_bool_arr = np.isin(prot_arr, combinations[i])
        vamp_struc_df_filter = vamp_struc_df[prot_bool_arr & selection_arr]
        calc_substitution_matrix(vamp_struc_df_filter,f"{outdir_pre}/{i}")
        
        
def prepare_data(wcn_cut, rasa_cut, use_monomers=False):
    """
    Read vampseq data and structural data for all proteins and make 
    selection_dict, which that for each structural feature contains 
    a boolean array specifying the indices of variants with that feature 
    in vamp_struc_df. 
    Select here whether to run for monomer or dimer data (which only changes 
    the structure features for NUDT15 and ASPA), and select cutoff values 
    for wcn and rasa for defining buried and exposed residues. 
    """
    # define proteins
    proteins = ['PTEN','TPMT','CYP2C9','NUDT15','ASPA','PRKN']

    # read vampseq data
    vampseq_df = pd.read_csv("../data/vampseq_data/processed/vampseq_soluble_df.csv")

    # read structure data (with NUDT15 and ASPA dimers) derived from AF structures
    if use_monomers == True:
        filename_add = "_monomers"
    else:
        filename_add = ""
    structure_df = pd.read_csv(f"../output/structure_features/structure_features{filename_add}_AF.csv") 
    structure_df = structure_df.rename(columns={"residue_number":"resi"})

    # merge vamp and structure dataframes
    vamp_struc_df = vampseq_df.merge(structure_df,how="left",on=["prot","resi"])

    # extract arrays with features from merged dataframe
    prot_arr = vamp_struc_df.prot.values 
    wcn_arr = vamp_struc_df.WCN.values
    rasa_arr = vamp_struc_df.rASA.values
    secondary_structure_arr = vamp_struc_df.secondary_structure.values

    # define structural features and create bool arrays 

    # buried and exposed
    buried_arr = (wcn_arr >= wcn_cut) & (rasa_arr <= rasa_cut)
    exposed_arr = (wcn_arr < wcn_cut) | (rasa_arr > rasa_cut)

    # helix, strand and loop
    helix_types = ["H","G","I"]
    strand_types = ["E","B"]
    loop_types = ["T","S","-"]
    alpha_helix = ["H"]

    helix_arr = np.isin(secondary_structure_arr, helix_types)
    strand_arr = np.isin(secondary_structure_arr, strand_types)
    loop_arr = np.isin(secondary_structure_arr, loop_types)
    alpha_helix_arr = np.isin(secondary_structure_arr, alpha_helix)

    # global is a special case where all should just evaluate to True
    global_arr = prot_arr.copy()
    global_arr[:] = True 

    # make dict with assays for each feature/environment
    selection_dict = dict(zip(["global","buried","exposed","helix","strand","loop","alpha_helix"],
                              [global_arr,buried_arr,exposed_arr,helix_arr,strand_arr,loop_arr,alpha_helix_arr]))
    
    return vamp_struc_df, selection_dict        


def scan_cutoff_values(wcn_cut_list, rasa_cut_list):
    """
    Calculate exposure-based substitution matrices for all c
    ombinations of wcn_cut and rasa_cut values in wcn_cut_list and
    rasa_cut_list. Use dimers structures for NUDT15 and ASPA. 
    """
    # loop over cutoff values
    for wcn_cut, rasa_cut in itertools.product(wcn_cut_list,rasa_cut_list):
    
        # prep matrix calc input data
        vamp_struc_df, selection_dict = prepare_data(wcn_cut, rasa_cut, use_monomers=False)

        # calculate matrices for individual structural features
        for feature in ["buried","exposed"]:

            # make output dirs if not already existing
        
            if not os.path.isdir(f"../output/substitution_matrices/dimers_gridsearch/{wcn_cut}_{rasa_cut}"):
                os.mkdir(f"../output/substitution_matrices/dimers_gridsearch/{wcn_cut}_{rasa_cut}") 
                
            if not os.path.isdir(f"../output/substitution_matrices/dimers_gridsearch/{wcn_cut}_{rasa_cut}/{feature}"):
                os.mkdir(f"../output/substitution_matrices/dimers_gridsearch/{wcn_cut}_{rasa_cut}/{feature}") 

            # calc matrices
            calc_matrices_for_selection(vamp_struc_df, selection_dict[feature],
             f"dimers_gridsearch/{wcn_cut}_{rasa_cut}/{feature}")


def calc_all_matrices(wcn_cut, rasa_cut, use_monomers=False): 
    """
    Calculate all substitution matrices with buried/exposed
    residue classifications based on wcn_cut and rasa_cut.
    """
    # define proteins
    proteins = ['PTEN','TPMT','CYP2C9','NUDT15','ASPA','PRKN']
    
    # prepare matrix calc input data
    vamp_struc_df, selection_dict = prepare_data(wcn_cut, rasa_cut, use_monomers=use_monomers)

    if use_monomers == True:
        struc_dir_name = "monomers"
    else:
        struc_dir_name = "dimers"
    
    # calculate matrices for individual structural features
    # (Global, Exposure, Secondary structure)
    for feature in selection_dict:

        # make output dir if not already existing
        if not os.path.isdir(f"../output/substitution_matrices/{struc_dir_name}/{feature}"):
            os.mkdir(f"../output/substitution_matrices/{struc_dir_name}/{feature}")  

        # calc matrices
        calc_matrices_for_selection(vamp_struc_df, selection_dict[feature], f"{struc_dir_name}/{feature}")

    # and then calculate matrices for combinations of structural features
    # (Combined)
    for feature_1,feature_2 in itertools.product(["buried","exposed"],["helix","strand","loop","alpha_helix"]):

        # make output dir if not already existing
        if not os.path.isdir(f"../output/substitution_matrices/{struc_dir_name}/{feature_1}_{feature_2}"):
            os.mkdir(f"../output/substitution_matrices/{struc_dir_name}/{feature_1}_{feature_2}")

        # calc matrices
        calc_matrices_for_selection(vamp_struc_df, selection_dict[feature_1] & selection_dict[feature_2], 
                                    f"{struc_dir_name}/{feature_1}_{feature_2}")

    # create dict with protein name combinations
    combinations = {}
    for i in [1,2,3,4,5]:
        count = 0
        for j in itertools.combinations(proteins,i):
            combinations[str(i)+'_'+str(count)] = list(j)
            count += 1

    # save dict with protein combinations        
    f = open(f'../output/substitution_matrices/{struc_dir_name}_combinations/combinations_dict.pkl', 'wb')
    pickle.dump(combinations, f)
    f.close()   

    # calculate matrices for all protein combinations for residues in exposed and buried environments
    for feature in ["buried","exposed"]:

        if not os.path.isdir(f"../output/substitution_matrices/{struc_dir_name}_combinations/{feature}"):
            os.mkdir(f"../output/substitution_matrices/{struc_dir_name}_combinations/{feature}")

        calc_matrices_for_selection_all_combinations(vamp_struc_df, selection_dict[feature], 
                                                     combinations, f"{struc_dir_name}_combinations/{feature}")
    
    # calculate matrices for all protein combinations for combined exposure and secondary structure environments
    for feature_1,feature_2 in itertools.product(["buried","exposed"],["helix","strand","loop"]):
    
        if not os.path.isdir(f"../output/substitution_matrices/{struc_dir_name}_combinations/{feature_1}_{feature_2}"):
            os.mkdir(f"../output/substitution_matrices/{struc_dir_name}_combinations/{feature_1}_{feature_2}")

        calc_matrices_for_selection_all_combinations(vamp_struc_df, selection_dict[feature_1] & selection_dict[feature_2], 
                                                     combinations, f"{struc_dir_name}_combinations/{feature_1}_{feature_2}")
    
        
if __name__ == "__main__":
    
    # run scan of wcn and rasa cutoff combinations
    # to find optimal values
    wcn_cut_list = [0,5,10,12.5,15,17.5,20]
    rasa_cut_list = [0,0.1,0.2,0.3,0.4,0.5]
    scan_cutoff_values(wcn_cut_list, rasa_cut_list)
    
    # run all substitution matrix calculations 
    # with optimal rasa and wcn cutoff values
    # (everything in all folders except dimers_gridsearch is based on the wcn_cut/rasa_cut given here)
    wcn_cut = 0
    rasa_cut = 0.1
    calc_all_matrices(wcn_cut, rasa_cut, use_monomers=False)
    calc_all_matrices(wcn_cut, rasa_cut, use_monomers=True)
