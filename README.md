# startmrca
This package implements a hidden Markov model to estimate the time at which a beneficial
allele begins increasing in frequency. More specifically, we use assume a "star" genealogy 
among carriers of the benefical allele to estimate time to the most recent common ancestor
(TMRCA). See [Smith et al](https://www.biorxiv.org/content/early/2016/08/24/071241) 
for details of the model and methods.

## Setup
In R, install the [latest version of the
startmrca package](https://github.com/jhavsmith/startmrca/releases) using these commands:
    
   ```R
   library(devtools)
   install_github("jhavsmith/startmrca@v0.6-1")
   library(startmrca)
   ```
