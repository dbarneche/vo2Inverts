# Do low oxygen environments facilitate marine invasions? Relative tolerance of native and invasive species to low oxygen conditions

This repository contains code and data needed to reproduce the article:

**Lagos ME, Barneche DR, White CR, Marshall DJ**, Do low oxygen environments facilitate marine invasions? Relative tolerance of native and invasive species to low oxygen conditions. *Global Change Biology*. doi:[10.1111/gcb.13668](http://onlinelibrary.wiley.com/doi/10.1111/gcb.13668/full) (accepted 2017-02-08) 

## Instructions

All analyses were done in `R`. To compile the paper, including figures and tables we use the [remake](https://github.com/richfitz/remake) package for R. You can install remake using the `devtools` package:

```r
devtools::install_github("richfitz/remake", dependencies=TRUE)
```
(run `install.packages("devtools")` to install devtools if needed.)

The `remake` package also depends on `storr`, install it like this:
```r
devtools::install_github("richfitz/storr", dependencies=TRUE)
```

You also need to install the software `JAGS`, as per these instructions [here](https://sourceforge.net/projects/mcmc-jags/?source=typ_redirect).

Next you need to open an R session with working directory set to the root of the project.

We use a number of packages, missing packages can be easily installed by remake:

```r
remake::install_missing_packages()
```

To generate all figures, analyses, and tables, simply do:

```r
remake::make()
```

All output will be automatically placed in a directory called `output` (it is going to be automatically created for you).

Also notice that the Bayesian analysis in this paper will take many hours to run on a regular computer.

If you find remake confusing and prefer to run plain R, you can use remake to build a script `build.R` that produces a given output, e.g.

```r
remake::make_script(filename="build.R")
```

### The paper can be reproduced using the following software and associated packages:
```
R version 3.3.2 (2016-10-31)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X El Capitan 10.11.6

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] R2jags_0.5-7      rjags_4-6         coda_0.18-1       plyr_1.8.4        LoLinR_0.0.0.9000

loaded via a namespace (and not attached):
[1] zoo_1.7-13       parallel_3.3.2   abind_1.4-5      Rcpp_0.12.9      grid_3.3.2       lmtest_0.9-34    boot_1.3-18      lattice_0.20-34  R2WinBUGS_2.1-21
```

### Please report if you run into problems or spot a bug on the code:
d13g0 DOT b4rn3ch3 AT m0n4sh DOT 3du (replace the 0 for o, 1 for i, 3 for e, 4 for a)  

### How to download this project for people not familiar with GitHub:  
* on the project main page on GitHub, click on the green button `clone or download` and then click on `Download ZIP`  
