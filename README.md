# Measuring temporal dynamics of metacommunities
 This repository contains data and code necessary to reproduce results shown in the paper:    
 ```
 Li HD, Holyoak M, Xiao ZS*. 2023. Disentangling spatiotemporal dynamics in metacommunities through a 
 species-patch network approach. Ecology Letters. (In revision)
```
Contact: *lihd AT ioz.ac.cn* (Hai-Dong Li)

## Metacommunity simulation
To simulate the metacommunity dynamics, we used the R package ["mcomsimr"](https://github.com/plthompson/mcomsimr) developled by:
```
Patrick L. Thompson, Laura Melissa Guzman, Luc De Meester, Zsófia Horváth, Robert Ptacnik, 
Bram Vanschoenwinkel, Duarte S. Viana, & Jon M. Chase. 2020. A process based framework for
metacommunity ecology. Ecology Letters. 23:1314-1329. (https://doi.org/10.1111/ele.13568)
```
Please note that some of the code in this folder was written to use parallel computing and run on a cluster. The simulation proccess was very very slow using R. One can use `Julia` for speeding up (see Thompson's EL paper: https://github.com/plthompson/Meta_com_framework).

## Case studies
The `Hoja mammal metacommunity data` can be accessed on Figshare: (https://doi.org/10.6084/m9.figshare.12978200.v1); the `Lindholm's vascular aquatic macrophyte metacommunity data` can be accessed Dryad Digital Repository: (https://doi.org/10.5061/dryad.t1g1jwsxv); the `Horvath's metacommunity data` can be accessed on Figshare
Repository: (https://doi.org/10.6084/m9.figshare.7823726).          
We provided R codes to analyze the species-patch networks case by case.

## Sensitivity analysis
The `Horvath's metacommunity data` can be accessed on Figshare Repository: https://doi.org/10.6084/m9.figshare.7823726.
We provided R codes to analyze sensitivity of both network structure of species-patch networks and beta diversity of speices-links to sampling efforts.
