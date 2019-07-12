## Multiple Mediation Analysis with Nonlinear Models

<img src="man/figs/dag_complex2.png" align="right" width="240" />


Fit mediation models to assess the causal impact of **x** on **y** mediated through **m**. Implements models and algorithms described in Long et al. 2019+.

### Installation

Install from github using R `devtools`

``` r
## easy access to vignettes, takes longer and requires suggests packages
devtools::install_github("longjp/mediateR", build_opts = c("--no-resave-data", "--no-manual"))

## faster, no vignettes
devtools::install_github("longjp/mediateR")
```

or in a terminal

``` r
git clone https://github.com/longjp/mediateR.git
R CMD INSTALL mediateR
```

### Examples

If `devtools::install_github` with `build_opts` was used, you can access vignettes with

``` r
browseVignettes("mediateR")
```

Otherwise download the `.html` files in the `\vignettes` folder, view in a browser.

### Citation and Contact

This code was developed for and described in Long et al. 2019+. Please cite this work if you use the templates for any publications. Email jplong@mdanderson.org with questions or bug reports.
