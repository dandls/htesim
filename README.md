# htesim - Simulation Study Framework for Heterogeneous Treatment Effect Estimation 

The study frameworks are based on the following papers: 

<!--- * Dandl S, Hothorn T, Heidi S (2020). Divide or Unite? A Comparison of Two Strategies to Random Forests-type Heterogeneous Treatment Effect Estimation. --->
* Wager S, and Athey S (2018). "Estimation and Inference of Heterogeneous Treatment Effects using Random Forests". Journal of the American Statistical Association, 113(523).
* Nie X, Wager S (2020). “Quasi-Oracle Estimation of Heterogeneous Treatment Ef-
fects.” Biometrika. ISSN 0006-3444. doi:10.1093/biomet/asaa076. Asaa076,
https://academic.oup.com/biomet/advance-article-pdf/doi/10.1093/biomet/asaa076/33788449/asaa076.pdf, 
URL https://doi.org/10.1093/biomet/asaa076.

## Project Status

Code is now on Github, our publication will be available soon. 

## Installation

You can install the github version, using `remotes`:

```r
remotes::install_github("susanne-207/htesim")
```

## Example 
The following code imitates Setup A of Nie and Wager (2020)  
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

MIT License
