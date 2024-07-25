# Three Open Questions in Polygenic Score Portability

Joyce Y. Wang<sup>1</sup>, Neeka Lin<sup>1</sup>, Michael Zietz<sup>2</sup>, Jason Mares<sup>3</sup>, Vagheesh M. Narasimhan<sup>1,4</sup>, Paul J. Rathouz<sup>4,5</sup> and Arbel Harpak<sup>1,5,+</sup>

<sub><sup>1</sup> Department of Integrative Biology, The University of Texas at Austin, Austin, TX</sub>

<sub><sup>2</sup> Department of Biomedical Informatics, Columbia University, New York, NY</sub>

<sub><sup>3</sup> Department of Neurology, Columbia University, New York, NY</sub>

<sub><sup>4</sup> Department of Statistics and Data Science, The University of Texas at Austin, Austin, TX</sub>

<sub><sup>5</sup> Department of Population Health, The University of Texas at Austin, Austin, TX</sub>

<sub><sup>+</sup> Correspondence should be addressed to A.H. (arbelharpak@utexas.edu)</sub>

## Abstract

A major obstacle hindering the broad adoption of polygenic scores (PGS) is their lack of “portability” to people that differ—in genetic ancestry, environmental exposures or other characteristics—from the GWAS samples in which genetic effects were estimated.  Here, we use the UK Biobank (UKB) to measure the change in individual-level PGS prediction accuracy as a continuous function of genetic dissimilarity to the GWAS sample ("genetic distance"). Our results highlight three major gaps in our understanding of PGS portability. First, prediction accuracy is extremely noisy at the individual level and not well predicted by genome-wide genetic distance. Second, mean trends of portability vary across traits. For several immune-related traits, prediction accuracy drops near zero quickly even at intermediate levels of genetic dissimilarity to the GWAS sample. This quick drop may reflect immune associations being more ancestry-specific than other trait associations. Third, we show that even qualitative trends of portability can depend on the measure of prediction accuracy used. For triglyceride count, a measure of prediction accuracy at the individual level (reduction in mean squared error) increases with genetic distance. Together, our results show that portability cannot be understood through global ancestry groupings alone. There are other, understudied drivers of portability, such as characteristics specific to the evolution of the trait and its genetic architecture, environmental determinants, and the construction of the polygenic score. Addressing these gaps can aid in the development and application of PGS and inform more equitable genomic research.

## Installation

We recommend creating a [conda](https://docs.conda.io/projects/conda/en/stable/) environment for running the code.

```
git clone https://github.com/harpak-lab/Portability_Questions
cd Portability_Questions
conda env create -f environment.yml
```

## Modification

Before execution, the directories contained in the scripts need to be modified so that they point to your directories.

## Execution

Execute all the bash scripts ending with `.sh` with `bash <script_name.sh>` in the following order:

1. `00_make_directories.sh`
2. `01a_extract_data_fields.sh`
3. `01b_filter_individuals_job.sh`
4. `01d_filter_genotype_files.sh`
5. `02_prepare_covariates_phenotypes.sh`
6. `03_gwas.sh`
7. `04a_clumping.sh`
8. `04e_after_clumping.sh`
9. `05a_pc_dist_fst.sh`
10. All the scripts created under `temp_fst_path`
11. `05e_find_best_num_pc.sh`
12. `05h_ukb_kgp_pca.sh`
13. `05j_pc_dist_fst_plots.sh`
14. `06_compute_prs.sh`
15. `07_group_ind_level_pred.sh`
16. `08a_prepare_close_far_pca.sh`
17. `08c_close_far_pca.sh`
18. `08d_prepare_close_far_gwas.sh`
19. `08f_gwas_close.sh` and `08f_gwas_far.sh`
20. `08h_calc_heterozygosity.sh`
21. `08k_compare_effect_sizes_heterozygosity_var_pgs_plots.sh`
