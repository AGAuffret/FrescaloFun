# FrescaloFun
This repository contains a .R file containing functions that will allow you to take an output folder from a Frescalo analysis (implemented using the `sparta` package), and use it to estimate recorder effort per grid cell, estimate and plot different (but quite similar) metrics of species' regional occupancy change, and calculate species' probability of occurrence in a given grid cell during a given time period. For instructions on running the functions, including example data, please see the Rmarkdown tutorial.

## Authors
[Alistair Auffret](mailto:alistair.auffret@slu.se) and Matilda Arnell

## Links and references
[Hill (2012)](https://doi.org/10.1111/j.2041-210X.2011.00146.x) - The original paper describing the Frescalo method.

[sparta](https://github.com/BiologicalRecordsCentre/sparta) - A common R package for running the Frescalo algorithm. Our functions require outputs from sparta.

[Pescott et al. 2022](https://doi.org/10.1016/j.ecolind.2022.109117) - Paper describing nice ways to estimate and visualise trend estimates and their uncertainty. Our functions for estimating and plotting trends include some of the methods described ehre.

[sacrevert](https://github.com/sacrevert) - Oli Pescott's github, invaluable in helping us actually understand frescalo.

[Fox et al. 2014](https://doi.org/10.1111/1365-2664.12256) - A paper looking at trends in British moths using Frescalo, and whose method for estimating trends is included in one of our functions.

[Eichenberg et al. 2021](https://doi.org/10.1111/gcb.15447) - A paper that maps probability of occurrence for species over time and produces a broad measure of species richness, which is included in one of our functions.
