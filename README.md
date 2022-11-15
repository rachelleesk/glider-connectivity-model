# glider-connectivity-model
[![DOI](https://zenodo.org/badge/565328138.svg)](https://zenodo.org/badge/latestdoi/565328138)
## About

This set of `rmd` files demonstrates the connectivity analysis for Sunda colugo in:
>  *Journal reference*

## Model

The current model supports the examination of connectivity between trees and vertical mitigation structures (e.g., glide poles). It also includes a genetic algorithm to determine optimal locations for the installation of vertical mitigation structures.

Information required for the model:
-	A map of trees with height and GPS locations
-	Digital elevation model
-	Glide ratio for studied species*

*In our publication, we used a glide ratio model instead of a single-value estimation. The glide ratio model is non-linear, and the glide ratio was found to increase with the horizontal glide distance.

## Data

The `data` set included in this repository is a mock data generated to demonstrate the workings of the code.

## Level of Support

The model was not built as an R package and we are simply sharing this with the community for wider application. While we do not provide any guarantees of support, the community is welcome to submit issues. We also welcome researchers and/or developers to take this further to the next level and develop it as a package for greater adaptability to other use cases. 

## Contact

Rachel S.K. Lee (rachelleesk@gmail.com)
