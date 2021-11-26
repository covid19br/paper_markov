# paper_markov

This repository implements the code used in the paper "Modelling optimal
vaccination strategies against Covid-19 in a context of Gamma variant
predominance in Brazil" by Ferreira et al.

[Preprint](https://www.medrxiv.org/content/10.1101/2021.11.19.21266590v1)

## Repository structure

-   functions/ contains (most) the functions used to prepare data,
    simulate and plot the results. The important files are:

    -    markov_all_vac.R: contains the model

    -    epi_params_all_vac.R: contains the epidemiologic parameters and
        initial structuring of the model

    -   vac_params_all_vac.R: contains the vaccination parameters of the
        model (but not really the effectiveness parameters)

    -   sensitivity_analysis_C1_2.R: contains the actual distribution of
        parameters used in the results

    -   SA\_{coverage, ideal, deaths, production}.R: contains the script
        to generate the results in the model. Notice that these scripts
        run several samples of the model and might take a while to run.
        To avoid this, these scripts provide the lines to load the
        tabulated results.

    -   run_all_vac.R: contains a simple script to plot a single result
        of the model

-   DATA/ contains the data used to structure the model, from age
    distribution to vaccine doses expected.

-   results/ contains the tabulated results.

-   plots/ contains the figures generated with the results

### Important dependencies

Tested in Ubuntu 64-bits 20.04. R version 4.1.2 (but expect to work with
version \> 3.6).

Packages Tidyverse, Matrix and doParallel are used.
