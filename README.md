# htesim - Simulation Study Framework for Heterogeneous Treatment Effect Estimation 

The study frameworks are based on the following papers: 

* Wager S and Athey S. Estimation and inference of heterogeneous treatment effects using random
forests. Journal of the American Statistical Association 2018. 113(523): 1228–1242.
Athey S, Tibshirani J and Wager S. Generalized random forests. The Annals of Statistics 2019.
47(2): 1148–1178.
* Nie X and Wager S. Quasi-oracle estimation of heterogeneous treatment effects. Biometrika 2021.
108: 299–319.

The package also includes the code for reproducing all findings of the paper in the respective folder [inst/empeval](https://github.com/dandls/htesim/tree/master/inst/empeval): 

Dandl S, Hothorn T, Seibold H, Sverdrup E., Wager S., Zeileis A. (2022). What Make Forest-based Heterogeneous Treatment
Estimators Work?. 


## Project Status

Code is now on Github, our publication will be available soon. 

## Installation

You can install the github version, using `remotes`:

```r
remotes::install_github("susanne-207/htesim")
```

## Example 
The following code immitates Setup A of Nie and Wager (2020)  
```r 
# Initialize manual for data generating process
# Functions for treatment effect, prognostic effect and treatment propensity 
dg <- dgp(t = tF_div_x1_x2, m = mF_sin_x1_x5, p = pF_sin_x3, model = "normal", xmodel = "unif")

# Simulate data 
sdg <- simulate(dg, nsim = 1000L, dim = 12) 
head(sdg) 
```

## Citation

* not yet determined 

## License

GPL-2 | GPL-3
