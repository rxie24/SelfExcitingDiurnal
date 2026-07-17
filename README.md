## Overview
This repository contains the code and summary data used to generate the tables and figures for the paper "Modeling Self-Exciting Diurnal Processes with Applications to Smartphone Use in Populations with Affective Disorders" by Ryan Xie, Sage Rush, Theodore D. Satterthwaite, and Ian Barnett.

Base functions needed to run both the simulation and real data analysis are in the "Base_Functions" folder. Summary data used in the real data analysis are stored in the "Data" folder. R scripts for fitting the Circadian Hawkes process and the Circadian Inverse Self-Exciting model to every individual's data are in the "CHP_Estimates" and "CISE_Estimates" folders respectively. Every subsequent folder houses the code used to generate their specific figure, table(s), or confidence intervals. To generate the results for figures 3 and 4, tables 4 and 5, and the parametric bootstrap confidence intervals, you first need to fit and save the respective model parameter estimates, in which code for doing so is in the "CHP_Estimates" and "CISE_Estimates" folders.

## Citation (APA)
Xie, R., Rush, S., Satterthwaite, T. D., & Barnett, I. (2026). Modeling self-exciting diurnal processes with applications to smartphone use in populations with affective disorders. Journal of Applied Statistics, 1-19.
