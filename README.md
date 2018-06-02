# offshore_mpc

This repository contains the source for the paper “Swiss Bank Savings versus Piggy Bank Savings: The Effect of Offshore Wealth on Estimates of Inequality and the Marginal Propensity to Consume”.

## Overview

Survey data containing wealth estimates at the household level is often used in conjunction with heterogeneous agent models to study the impact wealth inequality can have on an economy, including through the marginal propensity to consume (MPC). Unfortunately, these data often suffer from under-representation and differential non-response from wealthy households near the top of the distribution, leading to biased wealth distribution estimates. Using offshore wealth data from Alstadsæter, Johannesen, and Zucman (2018), we propose a novel method for mitigating these issues and apply our methodology to the ECB's Household Finance and Consumption Survey (HFCS). We re-calibrate a consumption-savings model proposed by Carroll, Slacalek, and Tokuoka (2014) using wealth distributions estimated with and without offshore. We find that correcting for offshore wealth leads to important changes in the predicted average MPC and distribution of MPCs for several countries. 

## File organisation

* **R/** - R source code
* **data/generated/** - data generated for or by calibrations
* **data/HFCS_UDB_1_3_ASCII/** - should be added for data summaries (see below)
* **deliverables/** - pre-compiled HTML files summarizing key results

## Requirements

* R 3.5
* RStudio 1.1 (recommended)
* R packages listed in **R/00_libraries.R**
* Python 3.5
* Python pakckage seaborn (only for MPC violin plot) 


## Summary of HFCS data

To estimate wealth distributions using HFCS data, access to the data needs to be requested from the [ECB](https://www.ecb.europa.eu/pub/economic-research/research-networks/html/researcher_hfcn.en.html).  Once authorization has been granted, the text-only source files should be placed in **data/HFCS_UDB_1_3_ASCII/**.  

To create summary statistics about wealth distributions estimated from the HFCS data, with and without the offshore adjustment, knit **data_summary.Rmd**.  

## Model calibration

To calibrate the model, change the parameters in and run **R/02_calibration.R**.  The country of interest (2 digit country code), target wealth distribution and input values can be changed at the beginng of the file.  A folder called **calibration_checkpoints/** should be created in the root directory to allow for progress to be saved.

## Summary of calibration results

The summary of the calibration results used in the paper can be created without having to re-calibrate the model by knitting **mpc_results.Rmd**.  The values for the Ginis and the solved models will be loaded from **ginis.RData** and **final_models.RData** (both inside **data/generated/**), respectively.