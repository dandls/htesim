This folder contains the code to run the benchmark study of Dandl et al. (2021). Experiments are based on batchtools. 

Files are encapsulated as following: Makefile -> run_study.R --> libs.R --> DGP.R --> def.R --> run.R
	* Makefile: Run `nohup make benchmark &` in your console to rerun the whole study. This will take several days depending on how many CORES are available (default CORES = 30L, setup in def.R). 
	* run_study.R: Code to run the study using the R package `batchtools`, called by the make command above.
	* def.R: Experimental setup to run study, called by run_study.R. Defines the number of trees, registry name, number of cores and name of results file.
	* create_plots_table.R: Code to reproduce plots and tables in Dandl et al. (2021). Calls helpers.R and setup.R. 
	* helpers.R: Helpers to create plots and tables. 
	* setup.R: Necessary packages & global plotting settings

