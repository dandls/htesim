This folder contains the code to run the simulation study of 
- Dandl et al. (2022): What Makes Forest-based Heterogeneous Treatment Effect Estimation Work?
--> paper1
- Paper 2 - Dandl et al. (2022): Heterogeneous Treatment Effect Estimation for Observational Data using Model-based Forests
--> paper2
 
For both papers, the experiments are encapsulated in the following files: run_study.R --> ../libs.R --> DGP.R --> def.R --> run.R
Experiments are based on the R package batchtools.
	* run_study.R: Code to run the study using the R package `batchtools`.
	* ../libs.R: Code to load necessary R packages.
	* DGP.R: Code to generate study setup (uses functions in htesim). 
	* def.R: Experimental setup to run study, called by run_study.R. Defines the number of trees, registry name, number of cores and name of results file.
	* run.R: Code to define a base model, to fit a forest, to get estimates and to return the MSE. Also includes the definition of batchtools methods. 

After the experiments have run, the following files allow for reproducing the figures and tables in both papers:  
	* create_plots_table.R: Code to reproduce plots and tables. Calls helpers.R and setup.R. 
	* helpers.R: Helpers to create plots and tables. 
	* setup.R: Necessary packages & global plotting settings
	

