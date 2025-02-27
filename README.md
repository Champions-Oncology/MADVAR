---
output: github_document
---


# MADVAR

<!-- badges: start -->
<!-- badges: end -->

A Variance Based Algorithm for Feature Selection in Biological data.

## Installation

You can install the development version of MADVAR like so:

```
library(devtools)  
install_github("Champions-Oncology/MADVAR")
```


## Examples  

This is an example for the function `madvar` for calculating a cutoff based on data variance and returning a filtered dataset accordingly:  
```
data(TCGAcolon_adenocarcinomaTPM)
filteredData <- madvar(data = TCGAcolon_adenocarcinomaTPM,
                       plot_density = F)
```
It is also possible to input a variance vector (using the `data` argument), for which a cutoff value is returned.  


You can also plot the cutoff over the distribution density in exploratory mode:
```
madvar(data = TCGAcolon_adenocarcinomaTPM,
       plot_density = T)

```

 
