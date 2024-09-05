import numpy as np
import pandas as pd
import mdtraj as md
import itertools
import pickle
from Bio.PDB.DSSP import make_dssp_dict


def make_dssp_df(protname,pdb_chain):
    """
    Read dssp output and return DataFrame with secondary 
    structure, ASA and RSA per residue. 
    """
    # read dssp output
    dssp_out_dict = make_dssp_dict(f'../output/dssp/{protname}.dssp')

    # get relevant info from dssp_dict
    dssp_dict = {}
    for res in list(dssp_out_dict[0].keys()):
        if res[0] == pdb_chain: # select only from specified chain
            dssp_dict[res[1][1]] = list(dssp_out_dict[0][res][0:3]) # get only secondary structure and ASA

    # create DataFrame with the selected DSSP output
    dssp_df = pd.DataFrame.from_dict(dssp_dict,orient='index',columns=['aa_structure','secondary_structure','ASA'])
    dssp_df = dssp_df.reset_index().rename(columns={'index':'residue_number'})

    # define maxASA (in Å^2) from Tien et al. (doi:10.1371/journal.pone.0080635)
    max_ASA = {'A':129.0, 'R':274.0, 'N':195.0, 'D':193.0, 'C':167.0, 'E':223.0, 'Q':225.0, 'G':104.0,
               'H':224.0, 'I':197.0, 'L':201.0, 'K':236.0, 'M':224.0, 'F':240.0, 'P':159.0, 'S':155.0,
               'T':172.0, 'W':285.0, 'Y':263.0, 'V':174.0}
    
    # normalise to relative solvent accessiblity
    maxASA_arr = np.array([max_ASA[aa] for aa in dssp_df.aa_structure.values])
    dssp_df['rASA'] = dssp_df.ASA.values / maxASA_arr
    
    return dssp_df


def WCN_chainA(protname, uniprot_seq):
    """
    Calculate WCN for every residue in chain 0, taking other chains 
    into account if they are present in the structure file. 
    Function works for monomers and multimers and always returns 
    the WCN for every residue in chain A.
    For non-glycine residues, inter-residue distances are evaluated 
    based on closest sidechain heavy atoms, and for pairs with glycine, 
    the closest distance between atoms are used. 
    """
    r0=7.0 # the radius for the calculation
    pdb=md.load(f'../data/pdb_files/crystal+alphafold_translated_structures/processed/{protname}.pdb') # load pdb with mdtraj
    topology=pdb.topology
    n_chains = topology.n_chains

    # define topologies with and without GLY 
    all_chains=topology.select('protein')
    all_chains_noGLY=topology.select(f'protein and not resname GLY') 

    # make arrays
    wcn=np.zeros((len(uniprot_seq)),dtype=float) 
    cm_adj=np.empty((len(uniprot_seq)*n_chains,len(uniprot_seq)*n_chains),dtype=float)
    cm_adj[:]=np.nan

    # calculate distances based first on closest atoms, then replace with sidechain-heavy where possible
    for pdb_select,scheme_e in zip([all_chains,all_chains_noGLY],['closest','sidechain-heavy']):

        # slice pdb
        pdb_slice=pdb.atom_slice(pdb_select)

        # make lists with residue indices and corresponding residue numbers
        residue_count = 0
        chain_count = 0
        residue_idx_list = []
        residue_number_list = []
        
        for chain in pdb_slice.top.chains:
            res_chain_list = []
            for residue in chain.residues:
                res_chain_list.append(residue_count)
                residue_count += 1
                residue_number_list.append(int(str(residue)[3:])+len(uniprot_seq)*chain_count)
            residue_idx_list.append(res_chain_list)
            chain_count += 1
        
        # make pair list
        pair_list = []
        for i in residue_idx_list[0]:
            for j in np.concatenate(residue_idx_list):
                if j > i + 2:
                    pair_list.append([i,j])

        # calculate distances between all residues 
        pdb_dist,pdb_rp=md.compute_contacts(pdb_slice,contacts=pair_list,scheme=scheme_e,periodic=False)
        cm=md.geometry.squareform(pdb_dist,pdb_rp)[0]

        # add calculated distances to cm_adj 
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                cm_adj[residue_number_list[i]-1,
                       residue_number_list[j]-1]=cm[i,j]

    # calculate wcn for all residues in chain A based
    # on distances to all residues in entire complex
    for i in range(len(uniprot_seq)):
        nan_flag=True
        for j in range(len(uniprot_seq)*n_chains):
            if np.isnan(cm_adj[i,j])!=True and cm_adj[i,j]!=0.0:
                nan_flag=False
                # calculate wcn given the distance between the two closest heavy atoms for each pair
                wcn[i]+=(1-(cm_adj[i,j]*10/r0)**6)/(1-(cm_adj[i,j]*10/r0)**12)
        if nan_flag==True:
            wcn[i]=np.nan

    return wcn


def make_structure_df(use_monomers=False):
    """
    Make DataFrame with structure all features. Use either only monomer 
    structures (use_monomers=True) or dimer structures for ASPA and NUDT15 
    and monomer structures for the rest (use_monomers=False).
    """
    # load dictionaries with uniprot names and sequences
    uniprotid_dict = pd.read_pickle("../data/uniprotid.pkl") 
    sequence_dict = pd.read_pickle("../data/sequences.pkl")
    
    # append _A to uniprot IDs for NUDT15 and ASPA if use_monomers=True
    if use_monomers==True:
        filename_add = "_monomers"
        for prot in ["NUDT15","ASPA"]:
            uniprotid_dict[prot] = uniprotid_dict[prot]+'_A'
    else:
        filename_add = ""
    
    structure_df_list = []

    for prot in ['PTEN','TPMT','CYP2C9','NUDT15','ASPA','PRKN']:

        print(prot)

        uniprot_id = uniprotid_dict[prot]
        sequence = sequence_dict[prot]

        # make dataframe with dssp output
        dssp_df = make_dssp_df(uniprot_id,"A")

        # calculate WCN for chain A
        wcn = WCN_chainA(uniprot_id, sequence)

        # add structural features to dataframe
        structure_df = pd.DataFrame()
        structure_df["residue_number"] = np.arange(1,len(sequence)+1,1)
        structure_df["WCN"] = wcn
        structure_df = structure_df.merge(dssp_df,how="left",left_on="residue_number",right_on="residue_number")
        structure_df["prot"] = prot

        # append
        structure_df_list.append(structure_df)

    structure_df_all = pd.concat(structure_df_list,ignore_index=True)
    structure_df_all.to_csv(f"../output/structure_features/structure_features{filename_add}_AF.csv",index=False)


if __name__ == "__main__":

    # make structure feature dataframe using NUDT and ASPA dimers
    # and monomers for all other proteins
    make_structure_df(use_monomers=False)
    
    # make structure feature dataframe using only monomers
    make_structure_df(use_monomers=True)