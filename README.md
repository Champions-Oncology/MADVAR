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
filteredData <- madvar(data = TCGAcolon_adenocarcinomaTPM,
                       plot_density = F)
```

You can also plot the cutoff in exploratory mode.
```
madvar(data = TCGAcolon_adenocarcinomaTPM,
       plot_density = T)

```

 
