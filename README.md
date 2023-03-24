# Measuring temporal dynamics of metacommunities
 This repository contains data and code necessary to reproduce results shown in the paper:    
 ```
 Li HD, Holyoak M, Xiao ZS*. 2023. Disentangling spatiotemporal dynamics in metacommunities through a 
 species-patch network approach. Ecology Letters. (In revision)
```
Contact: *lihd AT ioz.ac.cn* (Hai-Dong Li)

## Metacommunity simulation
We performed simulations to disentangle the relationships between both the structure of species-patch networks, and the temporal beta diversity of species-patch links and three core ecological processes. Two-stage simulations were used to teste if the temporal beta diversity of species-patch links partitioning helped to reveal the contribution of the processes that gains and losses of both species and patches to the metacommunity dynamics.    
Please note that some of the code in this folder was written to use parallel computing and run on a cluster. 

The simulation proccess was very very slow using R. One can use `Julia` for speeding up (see Thompson's EL paper: https://github.com/plthompson/Meta_com_framework).         

To simulate the metacommunity dynamics, we used the R package ["mcomsimr"](https://github.com/plthompson/mcomsimr) developled by:
```
Patrick L. Thompson, Laura Melissa Guzman, Luc De Meester, Zsófia Horváth, Robert Ptacnik, 
Bram Vanschoenwinkel, Duarte S. Viana, & Jon M. Chase. 2020. A process based framework for
metacommunity ecology. Ecology Letters. 23:1314-1329. (https://doi.org/10.1111/ele.13568)
```


## Case studies
We applied specie-patch network approaches to three long-term empirical metacommunity datasets, representing different scenarios. These included a mammal-patch system with biotic homogenization (`Arce-Pena et al., 2021`), an invertebrate species-pond system with habitat loss over time (`Horváth et al., 2019`), and a plant-lake system where environmental factors shaped species' spatial beta diversity over time (`Lindholm et al., 2020a`). We provided R codes to analyze the species-patch networks case by case.

The `Hoja mammal metacommunity data` (i.e., `Arce-Pena et al., 2021`) can be accessed on Figshare: (https://doi.org/10.6084/m9.figshare.12978200.v1); the `Lindholm's vascular aquatic macrophyte metacommunity data` can be accessed Dryad Digital Repository: (https://doi.org/10.5061/dryad.t1g1jwsxv); the `Horvath's metacommunity data` can be accessed on Figshare
Repository: (https://doi.org/10.6084/m9.figshare.7823726).          

## Sensitivity analysis
We used `Horvath's metacommunity data` to perform network sensitivity analysis. R codes to analyze sensitivity of both network structure of species-patch networks and beta diversity of speices-links to sampling efforts were provided.


The `Horvath's metacommunity data` can be accessed on Figshare Repository: https://doi.org/10.6084/m9.figshare.7823726.
