# _2024_Schulze_abundance-analysis

Data, scripts and notebooks for analysis in the manuscript:

To reproduce figures from manuscript, run the scripts in the following order to first calculate and plot average abundance score substitution matrices and then perform and plot predictions of variant abundance using the matrices:

- `/scripts/calc_structure_features.py` (calculate rASA and WCN for all residues in WT protein structures)

- `/scripts/calc_substitution_matrices.py` (calculate average abundance score substitution matrices for residues in different structural environments, using different combinations of variant abundance datasets and for different ways of classifying residues as exposed or buried)

- `/scripts/plot_substitution_matrices.py` 

- `/scripts/pred_from_matrices.py` (use average abundance score substitution matrices to predict variant abundance)

- `/scripts/plot_prediction_results.py` 

The plots shown in the manuscript that are not produced by running the above can be generated with the notebooks in `/notebooks`.

The `/data` directory contains protein sequences and UniProt IDs, the pooled VAMP-seq dataset used for meta-analysis, and structure files for all proteins. 

The structure files actually used in the analysis can be found in `/data/pdb_files/crystal+alphafold_translated_structures/processed/`.

The `/output/dssp` directory contains DSSP output for these structures, and features of these structures are stored in files in `/output/structure_features`. 

We have analysed VAMP-seq data from the following papers, and the files in `/data/vampseq_data/processed` thus contain data from these publications:

```

Matreyek KA, Starita LM, Stephany JJ, Martin B, Chiasson MA, Gray VE, Kircher M, Khechaduri A, Dines JN, Hause RJ, Bhatia S, Evans WE, Relling MV, Yang W, Shendure J, Fowler DM. Multiplex assessment of protein variant abundance by massively parallel sequencing. Nature Genetics. 2018 Jun; 50(6):874–882

Matreyek KA, Stephany JJ, Ahler E, Fowler DM. Integrating thousands of PTEN variant activity and abundance measurements reveals variant subgroups and new dominant negatives in cancers. Genome Medicine. 2021863 Dec; 13(1):165

Amorosi CJ, Chiasson MA, McDonald MG, Wong LH, Sitko KA, Boyle G, Kowalski JP, Rettie AE, Fowler DM, Dunham MJ. Massively parallel characterization of CYP2C9 variant enzyme activity and abundance. The American Journal of Human Genetics. 2021 Sep; 108(9):1735–1751.775

Suiter CC, Moriyama T, Matreyek KA, Yang W, Scaletti ER, Nishii R, Yang W, Hoshitsuki K, Singh M, Trehan A, Parish C, Smith C, Li L, Bhojwani D, Yuen LYP, Li Ck, Li Ch, Yang Yl, Walker GJ, Goodhand JR, et al. Massively parallel variant characterization identifies NUDT15 alleles associated with thiopurine toxicity. Proceedings of the National Academy of Sciences. 2020 Mar; 117(10):5394–5401

Grønbæk-Thygesen M, Voutsinos V, Johansson KE, Schulze TK, Cagiada M, Pedersen L, Clausen L, Nariya S, Powell RL, Stein A, Fowler DM, Lindorff-Larsen K, Hartmann-Petersen R. Deep mutational scanning reveals a correlation between degradation and toxicity of thousands of aspartoacylase variants. Nature Communications. 2024 May; 15(1):4026

Clausen L, Voutsinos V, Cagiada M, Johansson KE, Grønbæk-Thygesen M, Nariya S, Powell RL, Have MKN, Oestergaard VH, Stein A, Fowler DM, Lindorff-Larsen K, Hartmann-Petersen R. A mutational atlas for Parkin proteostasis. Nature Communications. 2024 Feb; 15(1):1541

```
