# Modeling virus epidemics taking into account the factor of social stress (COVID-19 SIR-ss model)

Materials for the global analysis of COVID epidemics SIRss model. 

### Introduction

The dynamics of epidemics hinge on the evolving behavior of individuals during an outbreak. In the initial stages, when awareness about the virus is low, people may not take preventive measures seriously. However, as alarms are raised and the outbreak intensifies, adherence to restrictions tends to increase, potentially curbing the spread of the epidemic. As time progresses, a phenomenon of fatigue and frustration with restrictions may set in, leading some individuals to discontinue compliance, especially if there's a decline in new cases. Following a period of respite, individuals may resume adhering to restrictions, yet this interval can expose them to the risk of a more potent second wave emerging during the pause. Notably, traditional SIR models often fail to predict the swift exit from the initial epidemic wave, emphasizing the importance of considering social dynamics. The emergence of a second wave is also intricately tied to social factors.

[Our published study](https://doi.org/10.1038/s41598-021-01317-z) delved into social stress through the lens of sociophysics. By amalgamating a dynamic SIR-type model with the classical triad of the general adaptation syndrome—alarm, resistance, exhaustion—we achieved a highly accurate depiction of the available statistical data.

In a subsequent study, a global analysis of statistical data encompassing 169 countries was conducted. Employing clustering and elastic maps for data dimensionality reduction, we scrutinized the distribution of parameters within the modified SIR model for these countries. This approach revealed a finite number of patterns characterizing the development of the COVID-19 pandemic. Armed with the understanding of these scenarios and their geographical correlations, local authorities can optimize their strategies for effectively combating epidemics until vaccines are developed.

### Global computation of parameters of SIRss COVID-19 epidemics model for 169 countries

Matlab code for the global computation of parameters of SIRss COVID-19 epidemics model for 169 countries can be found [here](https://github.com/lamhda/COVID_SIRss/blob/main/src/SIR_coeffs_to_2x169graphs.m). The resulting parameter table is [here](https://github.com/lamhda/COVID_SIRss/blob/main/data/results/Table_parameters.xlsx). A collection of plots in MATLAB fig format for individual countries
can be found [here](https://github.com/lamhda/COVID_SIRss/blob/main/data/results/2x169%20separate%20figures.zip) and [the MATLAB code for generating them](https://github.com/lamhda/COVID_SIRss/blob/main/src/SIR_coeffs_to_2x169separate_graphs.m).
The computations have been made based on [the datafile](https://github.com/lamhda/COVID_SIRss/blob/main/data/raw/owid-covid-data-10.05.2021.xlsx) obtained from [here](https://github.com/owid/covid-19-data/tree/master/public/data), accessed on 10.05.2021. 
The current version of the file can be downloaded [here](https://covid.ourworldindata.org/data/owid-covid-data.xlsx).

### Cluster analysis of the parameters of SIRss COVID-19 epidemics model for 169 countries

[This notebook](https://github.com/lamhda/COVID_SIRss/blob/main/COVID_geographical_var4.ipynb) contains the code for 
1. Visualizing individual model parameters on top of the geographical map
2. Visualizing combinations of parameters on top of the geographical map, using different color channels
3. Visualizing results of PCA analysis of parameter distributions on top of the geographical map, using different color channels
4. Visualizing results of elastic map method application on top of the geographical map
5. Performing clustering analysis of parameter distributions, using K-Means, DBScan, hierarchical ans spectral clustering, and applying clustering quality criteria
6. Visualizing clustering results on top of the geographical map

### References

Kastalskiy, I.A., Pankratova, E.V., Mirkes, E.M. et al. Social stress drives the multi-wave dynamics of COVID-19 outbreaks. Sci Rep 11, 22497 (2021). [https://doi.org/10.1038/s41598-021-01317-z]
