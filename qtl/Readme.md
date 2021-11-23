These are steps for running a simple eQTL analysis on the bulk RNA data where each 'visit' is run separately per AMP-PD cohort; ie simple style.

Create a JupyterLab instance using the Google Cloud Platform (GCP) AI Platform Notebooks interface within the PPMI WB RNA GCP project. For running the simple eQTL analysis with tensorQTL below I used a n1-standard-16 (16 vCPUs, 60 GB RAM) instance with 1024 GB data disk and 1 NVIDIA Tesla T4.

The GCP cloud console interface for the PPMI WB RNA project is here [AI Platform Notebooks](https://console.cloud.google.com/ai-platform/notebooks/list/instances?project=ppmi-wb-rna). The GCP documentation for here [AI Platform Notebooks documentation](https://cloud.google.com/ai-platform/notebooks/docs).

Once the instance has been created you should see the OPEN JUPYTERLAB link to access the JupyterLab instance.
Note: this JupyterLab is missing some handy extensions (such as TOC) and difficult to add them in this instance

Within JupyterLab open a Terminal
1. setup and data pulls
    - make a directory for the notebooks

        mkdir notebooks

    - pull the setup notebook

        gsutil cp gs://ppmi-wb-analysis/eqtl/notebooks/setup_and_data_pull.ipynb notebooks/

2. open the notebooks/setup_and_data_pull.ipynb and run, this will pull down the rest of the notebooks and input data such as genotypes, expression, gene annotations, and some sample info

3. run notebooks/foundin_format_expression_covariates.ipynb, which will pull the various known subject and sample covariates, clean them up as needed and save into a single covariate file, appropriate for this analysis, currently using:
     - AMP-PD subject info: amppd_demographicsPlus_2019_v1release_1015.csv
     - FOUNDIN-PD cell line info: cell_metadata.csv
     - Genetics PCs generated from PPMI WGS genotypes: foundin.freeze9.pca.eigenvec
     - estimated cell fractions from SCADEN deconvolution: rnab_cell_fracs_scaden.csv
     - rnaseq QC metrics: foundin_rnab_seqqc_metrics.txt
     - FOUNDIN-PD subject info: Expanded_overview_of_included_PPMI_samples_overview.csv
     - FOUNDIN-PD PD risk GRS: Expanded_overview_of_included_PPMI_samples_GRS.csv

4. run the notebooks/foundin_split_expression_by_day.ipynb, which will split the expression data (here using vst, this can be changed notebook variables) by days present in the data, this notebook has temp fix on renaming the few 'da' samples to 'd'

5. run notebooks/foundin_frmt_tensorqt_genos.ipynb, tensorqtl uses panda_plink to load the genotypes which basically reads plink bfiles, since we typically just use per chrom vcfs or plink pfiles need to convert to plink pfiles, and since cohort is pretty small and tensorQTL pretty fast go ahead and merge the per chrom pfiles into autosome bfile
    
- note since basically repeating prep and analysis over three separate timepoints individually I used [Papermill](https://papermill.readthedocs.io/en/latest/) automate the notebook runs instead of manually creating three separate experiment notebooks

6. Prep the expression data by running notebooks/foundin_expr_prep_runner.ipynb, this notebook uses papermill to run notebooks/foundin_day_expression_prep.ipynb where I just change to day variable (parameter), it will then create an output notebook for that run (experiment). This notebook, performs sex check from expression data, subsets the expression features based on detection rate (here 75%), quantile normalizes the data, and adjusts the expression features by coveriates generated from the expression features using UMAP projects from expression features excluding features from the bottom quartile of variance. The notebook formats the normalized and adjusted expression features matrix as a Hierarchical Data Format version 5 (HDF5) files and then also per chromosome plink2 formatted pheno files. 

    - note since in the above prep step since I generated and adjusted the expression data based on covariate generated from the data, this is one reason why I did the individual differention days separately as I expect this would be the primary covariate generated from the data together.

    - Review the expression data using the Papermill generated notebooks. For instance the notebook generated and run for the d25 expression data would be found notebooks/pm_gend_nbs/foundin_d25_expression_prep.ipynb. Here everthing should be fine but prior to running the eQTL analysis would want to verify how things looked as far as how many genes detected, all samples are expected sex, normalization doesn't look off, and how the adjustment for generated covariates changed the data (also includes some information on what known covariates, as they become available(?), are correlated with the generated ones). I often do look to see if the generated covaritates are going to remove signal related to variables that may be used in other analyses for future reference; ie disease status. Here for simple eQTL where detecting the effect of genotype on expression it is OK if covariate removes that signal so we are hopefully looking at closer to pure genotype effect. Of course this should be typically OK from eQTL at risk variances but would not hold for common highly penetrant variants, which we don't really have for neuroD.

7. Run the simple eQTL analysis, again using Papermill to run the days separately from same template notebook, by running the notebooks/foundin_cis_eqtl_runner.ipynb notebook. The primary outputs here the generated notebooks, such as notebooks/pm_gend_nbs/foundin_d65_cis_eqtl.ipynb, all analysis results per chromosome as hdf5 files, and all significant results based on B&H FDR at 0.05 again as hdf5 files. The actual eQTL analysis is run using Plink2 glm.

8. Also added a peek runner that looks at CHURC1, LRRK2, GBA, and SNCA eQTL results and then any eQTL that included PD Meta5 variants. CHURC1 is included as a sanity check and give an idea of context; typically if you run eQTL you should see this eQTL and it will probably be one of the largest you see. Local quick Manhattan's are generated for those four genes, where only a list of gene names is spill for the meta5 variants.

    - Note results on first peek are not super exciting, but this simple analysis is based on running FDR for full transcriptome analysis, next I'll run focused on DX risk variants only, then later for genes of interest (ie diff expression genes or eQTS genes). Basically p-value hacking to reduce tests for small sample size, so focus on stuff of interest.

9. To button everything up run the notebooks/finish_and_push_back.ipynb notebook to push results, info files, and generated notebooks back to the Google cloud storage (GCS) buckets

10. If completely done you should delete the compute instance from the GCP AI Platform Notebooks interface. Note if you feel like you need to do some more exploratory analysis or poking around you can just stop the instance instead of deleting and then come back later a start the instance again and jump back in. The instance disks are persistent, obviously memmory state isn't presevered, but data and tooling on the disk are. However while the cost is much reduced by not having the compute resources running we are still paying small fee in storage for that persistent disk and machine image so don't leave in place if not needed any longer. By the way sometimes helpful, you can stop the instance and change the machine type and restart so if need more or less cpu, memory, even gpu can do.


