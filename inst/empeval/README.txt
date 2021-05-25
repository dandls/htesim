This folder contains the code to run the benchmark study of Dandl et al. (2021). Experiments are based on batchtools. 

Files are encapsulated as following: Makefile -> run_study.R --> libs.R --> DGP.R --> def.R --> run.R
	* Makefile: Run `nohup make benchmark &` in your console to rerun the whole study. This will take up to 7 days depending on how many CORES are available (default CORES = 30L). 
	* run_study.R: Code to run the study using the R package `batchtools`, called by the make command above.
	* def.R: Experimental setup to run study, called by run_study.R.  
	* def_test.R: Experimental setup for a test run of the study, differs to def.R in its definition of number of trees, registry name, number of cores and name of resultsfile. It is called by run_study.R under `TEST = TRUE`. 
	* reduce_results.R: Code to preprocess results from batchtools registry. 
