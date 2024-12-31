# Three Open Questions in Polygenic Score Portability

Joyce Y. Wang<sup>1</sup>, Neeka Lin<sup>1</sup>, Michael Zietz<sup>2</sup>, Jason Mares<sup>3</sup>, Vagheesh M. Narasimhan<sup>1,4</sup>, Paul J. Rathouz<sup>4,5</sup> and Arbel Harpak<sup>1,5,+</sup>

<sub><sup>1</sup> Department of Integrative Biology, The University of Texas at Austin, Austin, TX</sub>

<sub><sup>2</sup> Department of Biomedical Informatics, Columbia University, New York, NY</sub>

<sub><sup>3</sup> Department of Neurology, Columbia University, New York, NY</sub>

<sub><sup>4</sup> Department of Statistics and Data Science, The University of Texas at Austin, Austin, TX</sub>

<sub><sup>5</sup> Department of Population Health, The University of Texas at Austin, Austin, TX</sub>

<sub><sup>+</sup> Correspondence should be addressed to A.H. (arbelharpak@utexas.edu)</sub>

Provided below are instructions and details for scripts used to generate the results and figures in ["Three Open Questions in Polygenic Score Portability"](https://www.biorxiv.org/content/10.1101/2024.08.20.608703v1).

## Installation

Install the software:

1. [<i>Plink 1.9</i> (Purcell, S. & Chang, C)](https://www.cog-genomics.org/plink/)
2. [<i>Plink 2.0</i> (Purcell, S. & Chang, C)](https://www.cog-genomics.org/plink/2.0/)

Download the [UK Biobank](https://www.ukbiobank.ac.uk/) (UKB) dataset, following their guidelines. The scripts also use the [1000 Genomes phase 3 dataset provided by Plink](https://www.cog-genomics.org/plink/2.0/resources), but it is not necessary to download it beforehand, as `05h_ukb_kgp_pca.sh` contains scripts for downloading it.

For running the scripts, we recommend creating a [conda](https://docs.conda.io/projects/conda/en/stable/) environment.

```
git clone https://github.com/harpak-lab/Portability_Questions
cd Portability_Questions
conda env create -f environment.yml
```

## Modification

Before execution, the directories contained in the scripts need to be modified so that they point to your directories.

## Execution

Execute the bash scripts ending with .sh with bash <script_name.sh>. Please see the details for each script in the following sections.

### Preparing the data

Execute the following files to filter and prepare the data:

1. `00_make_directories.sh`
2. `01a_extract_data_fields.sh` (make sure to edit the file so it's pointing to the correct UKB basket file)
3. `01b_filter_individuals_job.sh`
4. `01d_filter_genotype_files.sh`
5. `02_prepare_covariates_phenotypes.sh`

### GWAS

In the selection of the GWAS sample, we used the White British classification as provided by the UKB.

Execute the following files to perform GWAS, clumping, and thresholding:

1. `03_gwas.sh`
2. `04a_clumping.sh`
3. `04e_after_clumping.sh`

### Genetic distance calculations

The fixation index (<i>F<sub>st</sub></i>) is a natural metric, a single number, to measure the divergence between two sets of chromosomes and we considered using it to measure the distance between the pair of chromosomes of an individual and chromosomes in the GWAS sample. However, calculating Fst was computationally costly, so we used Euclidean distance in the PC space as a single number proxying genetic distance from the GWAS sample.

Execute the following files to calculate <i>F<sub>st</sub></i>:

1. `05a_pc_dist_fst.sh`
2. All the scripts created under `temp_fst_path`

Then, execute the following files to calculate Euclidean distance:

1. `05e_find_best_num_pc.sh` (creates <b>Fig. S1</b>)
2. `05h_ukb_kgp_pca.sh` (downloads [1000 Genomes phase 3 dataset provided by Plink](https://www.cog-genomics.org/plink/2.0/resources))
3. `05j_pc_dist_fst_plots.sh` (creates <b>Fig. 1</b>)

### PGS and evaluating PGS prediction accuracy

Execute the following file to calculate PGS:

1. `06_compute_prs.sh`

We evaluated PGS prediction accuracy at both the group level and individual level:

1. `07_group_ind_level_pred.sh` (creates <b>Fig. 2, S2-13</b>)

We compared the variance in squared prediction error explained for 8 raw measures: genetic distance, Townsend Deprivation Index, average yearly total household income before tax, educational attainment, which we converted into years of education, minor allele counts for SNPs with different magnitudes of effects (three equally-sized bins of small, medium, and large squared effect sizes, see Fig. S23), and minor allele counts of all SNPs:

1. `08a_prepare_for_ma_counts.sh`
2. `08b_calc_ma_counts.sh`
3. `08d_ind_pred_plots.sh` (creates <b>Fig. 3, S14-21</b>)

### Additional analyses on lymphocyte count

To understand why immunity-related traits like lymphocyte count have group-level prediction accuracy that drops near zero even at a short genetic distance, we performed additional analyses.

We first performed two additional GWASs and compared the allelic effects across the three GWASs:

1. `09a_prepare_close_far_pca.sh`
2. `09c_close_far_pca.sh`
3. `09d_prepare_close_far_gwas.sh`
4. `09f_gwas_close.sh` and `09f_gwas_far.sh`

We calculated heterozygosity at index SNPs as a function of genetic distance:

1. `09g_calc_heterozygosity.sh`

We examined the variance of PGS as a function of genetic distance:

1. `09j_compare_effect_sizes_heterozygosity_var_pgs_plots.sh` (creates <b>Fig. 4, S22</b>)

We estimated the heritability associated with each index SNP:

1. `10a_compare_heritability.sh` (creates <b>Fig. S23-24</b>)

### Portability of disease traits

We estimated the PGS portability of 3 disease traits at the group level.

For these disease traits, we ran GWAS, clumped and thresholded the SNPs, and calculated PGS:

1. `11a_gwas_disease.sh`
2. `11c_clumping_disease.sh`
3. `11f_after_clumping_disease.sh`
4. `11g_compute_pgs_disease.sh`

Then we estimated group level portability of the disease traits:

1. `11h_group_level_pred_disease.sh` (creates <b>Fig. 5</b>)

### Distribution of important variables in the dataset

We first plotted the distribution of Townsend deprivation index, household income, sex, age, and country as a function of genetic distance:

1. `12_townsend_income_sex_age_country.sh` (creates <b>Fig. S25-29</b>)

Then we plotted the correlation between <i>F<sub>st</sub></i> and PCs 1, 2, 3, and 40:

1. `13_pcs_vs_fst.sh` (creates <b>Fig. S30-33</b>)
