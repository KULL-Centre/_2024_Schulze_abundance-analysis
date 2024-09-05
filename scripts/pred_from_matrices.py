import numpy as np
import pandas as pd
import itertools
import pickle
import os


def classify_target_prot_residues(pred_prot, wcn_cut, rasa_cut):
    """
    Create dictionary that for each structural feature contains a 
    boolean array specifying for each residue if this residue belongs 
    to that specific feature class or not. 

    pred_prot : name of protein for which to perform predictions
    """
    # get structure data
    structure_df = pd.read_csv("../output/structure_features/structure_features_AF.csv") # run with NUDT15 and ASPA dimer structures
    structure_df = structure_df.rename(columns={"residue_number":"resi"})
    structure_df_prot = structure_df.groupby("prot").get_group(pred_prot) 
    
    # make array with residue numbers 
    resi_arr = structure_df_prot.resi.values

    # make arrays with residue structure features
    rasa_arr = structure_df_prot.rASA.values
    wcn_arr = structure_df_prot.WCN.values
    secondary_structure_arr = structure_df_prot.secondary_structure.values

    # define structural features and create bool arrays
    global_arr = np.zeros(len(resi_arr),dtype='bool')
    global_arr[:] = True
    
    # buried/exposed
    buried_arr = (wcn_arr >= wcn_cut) & (rasa_arr <= rasa_cut)
    exposed_arr = (wcn_arr < wcn_cut) | (rasa_arr > rasa_cut) 

    # helix, strand and loop
    helix_types = ["H","G","I"]
    strand_types = ["E","B"]
    loop_types = ["T","S","-"]
    
    helix_arr = np.isin(secondary_structure_arr, helix_types)
    strand_arr = np.isin(secondary_structure_arr, strand_types)
    loop_arr = np.isin(secondary_structure_arr, loop_types)

    # make dict with arrays for each feature/environment
    residue_class_dict = dict(zip(["global","buried","exposed","helix","strand","loop"],
                              [global_arr,buried_arr,exposed_arr,helix_arr,strand_arr,loop_arr]))
    
    for class1,class2 in itertools.product(["buried","exposed"],["helix","strand","loop"]):
        residue_class_dict[f"{class1}_{class2}"] = residue_class_dict[class1] & residue_class_dict[class2]
    
    return resi_arr, residue_class_dict


def create_matrix_dict(matrix_name, matrix_type): 
    """
    Make dictionary with substitution matrices for all structural features.

    matrix_name : either all, name of single protein, or leave out single protein.
    matrix_type : specify whether to use mean or median abundance score.
    """
    # make dict
    matrix_dict = {}

    # loop over single structure features and add matrices to dict
    for feature in ["global","buried","exposed","helix","strand","loop"]:

        matrix = pd.read_pickle(f"../output/substitution_matrices/dimers/{feature}/{matrix_name}/{matrix_type}.pkl")
        matrix_dict[feature] = matrix

    # loop over combined structure features and add matrices to dict
    for feature_1,feature_2 in itertools.product(["buried","exposed"],["helix","strand","loop"]):

        feature = feature_1+"_"+feature_2

        matrix = pd.read_pickle(f"../output/substitution_matrices/dimers/{feature}/{matrix_name}/{matrix_type}.pkl")
        matrix_dict[feature] = matrix
        
    return matrix_dict


def pred_from_matrix(resi_arr,resi_class_arr,matrix_dict,sequence):
    """
    Predict abundance scores for all possible single residue substitution variants of a 
    sequence using substitution matrices with average abundance score. Predictions are made
    based on the structural environment of each residue, given in resi_class_arr, by 
    looking up abundance score averages in the structure-context specific substitution matrices.

    resi_arr : array with residue numbers
    resi_class_arr : array specifying the structural class/environment for all residues
    matrix_dict : dictionary with substiution matrices for the different structural environments
    sequence : sequence of the protein for which predictions are performed
    """
    # make empty lists
    variants = []
    predictions = []

    # loop over residues
    for resi in resi_arr: 

        # determine reference amino acid type for resi
        aa_ref = sequence[resi - 1]

        # determine structural class/type for resi
        struc_type = resi_class_arr[resi - 1]              

        # perform prediction if struc_type is in matrix_dict
        if struc_type in list(matrix_dict.keys()):
        
            # look up predictions in matrix based on resi amino acid type and structural class
            pred_arr = matrix_dict[struc_type].loc[aa_ref].values

            # make array with variant amino acid types
            aa_var = matrix_dict[struc_type].loc[aa_ref].index

            # concatenate aa_ref, resi and aa_var data
            var_arr = (pd.Series([aa_ref]*20) + pd.Series([str(resi)]*20) + pd.Series(aa_var)).values

            assert len(pred_arr) == len(var_arr)

            # append
            variants.append(var_arr)                                
            predictions.append(pred_arr)

        # otherwise skip residue
        else:
            continue
            
    # concatenate data
    variants = np.concatenate(variants)
    predictions = np.concatenate(predictions)    
        
    return variants, predictions


def run_pred_dimers(pred_prot, matrix_name, matrix_type, wcn_cut, rasa_cut): 
    """
    Perform substitution matrix-based predictions for pred_prot using matrices 
    with name matrix_name (usually leave_out_pred_prot). 

    pred_prot : name of protein for which to perform predictions
    matrix_name : either all, name of single protein, or leave out single protein.
    matrix_type : specify whether to use mean or median abundance score.
    """
    # get sequence data
    sequence_dict = pd.read_pickle("../data/sequences.pkl")
    sequence = sequence_dict[pred_prot]

    # create dict with all matrices 
    matrix_dict = create_matrix_dict(matrix_name, matrix_type)

    # create dict with structure feature classifications for pred_prot
    resi_arr, residue_class_dict = classify_target_prot_residues(pred_prot, wcn_cut, rasa_cut) 

    # define structure types
    struc_type_dict = {"global":["global"],
                       "exposure":["buried","exposed"],
                       "secondary_structure":["helix","strand","loop"],
                       "combined":["buried_helix", "buried_strand", "buried_loop", 
                                   "exposed_helix", "exposed_strand", "exposed_loop"]}
    
    # make empty list to store predictions for each structure type
    pred_df_list = []
    
    # loop through structure types and do pred for each
    for struc_type in ["global","exposure","secondary_structure","combined"]: 
    
        # make array with residue classifications for structure type
        resi_class_arr = np.zeros(len(sequence),dtype="object")
        resi_class_arr[:] = np.nan
        for feature in struc_type_dict[struc_type]:
            resi_class_arr[residue_class_dict[feature]] = feature               
            
        # perform predictions   
        variants, predictions = pred_from_matrix(resi_arr, resi_class_arr, matrix_dict, sequence)

        # add predictions to dataframe
        temp_pred_df = pd.DataFrame()
        temp_pred_df["variant"] = variants
        temp_pred_df[struc_type] = predictions
        pred_df_list.append(temp_pred_df)
    
    # merge all dataframes in pred_df_list to have all predictions in single dataframe
    pred_df = pred_df_list[0]
    for i in range(1,len(pred_df_list)):
        pred_df = pred_df.merge(pred_df_list[i],how="outer",on=["variant"])    
    
    # save results
    pred_df.to_csv("../output/predictions_from_matrices/dimers/"+
                   f"pred_{pred_prot}_from_{matrix_name}_{matrix_type}.csv",index=False)

    
def run_pred_dimers_gridsearch(pred_prot, matrix_name, matrix_type):
    """
    Perform substitution matrix-based predictions for pred_prot
    using exposure-based matrices calculated with a range of
    wcn_cut and rasa_cut values to defined buried and exposed residues. 
    """
    # get sequence data
    sequence_dict = pd.read_pickle("../data/sequences.pkl")
    sequence = sequence_dict[pred_prot]
    
    # define cutoff values to try
    wcn_cut_list = [0,5,10,12.5,15,17.5,20]
    rasa_cut_list = [0,0.1,0.2,0.3,0.4,0.5]
   
    # initialise
    variants_list = []
    pred_df = pd.DataFrame()
    
    # loop over all wcn_cut and rasa_cut combinations and perform predictions
    for wcn_cut, rasa_cut in itertools.product(wcn_cut_list,rasa_cut_list):
    
        # make matrix dict for wcn_cut/rasa_cut combination
        matrix_dict = {}
        for feature in ["buried","exposed"]:
            matrix = pd.read_pickle("../output/substitution_matrices/dimers_gridsearch/"+
                                    f"{wcn_cut}_{rasa_cut}/{feature}/{matrix_name}/{matrix_type}.pkl")
            matrix_dict[feature] = matrix

        # classify residues according to wcn_cut and rasa_cut
        resi_arr, residue_class_dict = classify_target_prot_residues(pred_prot, wcn_cut, rasa_cut) 

        # make array with residue classifications based on residue_class_dict
        resi_class_arr = np.zeros(len(sequence),dtype="object") 
        resi_class_arr[:] = np.nan
        for feature in ["buried","exposed"]:
            resi_class_arr[residue_class_dict[feature]] = feature                   

        # perform prediction
        variants, predictions = pred_from_matrix(resi_arr, resi_class_arr, matrix_dict, sequence)
        
        variants_list.append(variants)
        
        pred_df[f"{wcn_cut}_{rasa_cut}"] = predictions

    # check that all var arrays are equal and add var column    
    for i in range(len(variants_list)):
        assert np.all(variants_list[0] == variants_list[i])
    pred_df.insert(loc=0, column="variant", value=variants_list[0])
    
    # save results
    pred_df.to_csv("../output/predictions_from_matrices/dimers_gridsearch/"+
                   f"pred_{pred_prot}_from_{matrix_name}_{matrix_type}.csv",index=False)


def run_pred_dimers_combinations(pred_prot, matrix_type, wcn_cut, rasa_cut): 
    """
    Perform substitution matrix-based predictions for pred_prot
    using exposure-based matrices or combined matrices. 
    Predictions are performed for all combinations of 
    proteins/datasets in which pred_prot does not occur.
    """
    # get sequence data
    sequence_dict = pd.read_pickle("../data/sequences.pkl")
    sequence = sequence_dict[pred_prot]

    # create dict with feature classifications for pred_prot
    resi_arr, residue_class_dict = classify_target_prot_residues(pred_prot, wcn_cut, rasa_cut)

    # define all structure types
    struc_type_dict = {"exposure":["buried","exposed"],
                       "combined":["buried_helix", "buried_strand", "buried_loop", 
                                   "exposed_helix", "exposed_strand", "exposed_loop"]}
    
    for struc_type in ["exposure","combined"]: 
    
        # define empty array with residue classifications
        resi_class_arr = np.zeros(len(sequence),dtype="object")
        resi_class_arr[:] = np.nan

        # add classifications to resi_class_arr
        for feature in struc_type_dict[struc_type]:
            resi_class_arr[residue_class_dict[feature]] = feature 

        # get combinations dict keys corresponding to matrices that pred_prot was not used to calculate
        combinations = pd.read_pickle("../output/substitution_matrices/dimers_combinations/combinations_dict.pkl")
        pred_combinations = []
        for i in combinations:
            if pred_prot not in combinations[i]: # perform pred using all combinations that pred_prot is not in
                pred_combinations.append(i)

        # define dataframe to store prediction results        
        pred_df = pd.DataFrame()

        # define variant list
        variants_list = []

        # loop through combinations without pred_prot and perform predictions
        for i in pred_combinations:

            # make matrix dict for combination i
            matrix_dict = {}
            for feature in struc_type_dict[struc_type]:
                matrix_dict[feature] = pd.read_pickle("../output/substitution_matrices/dimers_combinations/"+
                                                      f"{feature}/{i}/{matrix_type}.pkl")

            # perform prediction and add to dataframe
            variants, predictions = pred_from_matrix(resi_arr, resi_class_arr, matrix_dict, sequence)
            pred_df[i] = predictions

            # append
            variants_list.append(variants)

        # add variant column to dataframe    
        for i in range(len(variants_list)):
            assert np.all(variants_list[0] == variants_list[i])    
        pred_df.insert(loc=0, column="variant", value=variants_list[0]) 

        # save dataframe
        pred_df.to_csv("../output/predictions_from_matrices/dimers_combinations/"+
                       f"pred_{pred_prot}_{struc_type}_combinations_{matrix_type}.csv")    
      
    
if __name__ == "__main__":
    
    wcn_cut = 0
    rasa_cut = 0.1
    
    for matrix_type in ["mean","median"]:
    
        for pred_prot in ["PTEN","TPMT","CYP2C9","NUDT15","ASPA","PRKN"]:

            print(pred_prot)

            # perform all predictions using matrices calculated
            # by leaving out abundance data from pred_prot
            matrix_name = f"leave_out_{pred_prot}"

            # run all pred functions
            run_pred_dimers(pred_prot, matrix_name, matrix_type, wcn_cut, rasa_cut) 
            run_pred_dimers_combinations(pred_prot, matrix_type, wcn_cut, rasa_cut)
            run_pred_dimers_gridsearch(pred_prot, matrix_name, matrix_type)
