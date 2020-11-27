# SIR Modeling

<p align="center"><img src="https://www.lewuathe.com/assets/img/posts/2020-03-11-covid-19-dynamics-with-sir-model/sir.png"/></p>

`SIR_Modeling` is a project focuses on developing and examining the effects of different parameters on Susceptible-Infected-Recovered (SIR) models. The COVID-19 pandemic has led to the use of numerous epidemic models based upon SIR. Our goal is to understand the complexities of infectious diseases via 2 models: with and without intervention measures, as well as to underpin such models with real world data to understand their predictive capabilities and limitations.Home.

The models developed in R by [Linh Tang](https://linhtang.me/), [Britney He](https://github.com/britneyhe), and [Bowen Mince](https://github.com/minceb) under the mentorship of Professor Shona Kuiper as the final project for `STA310: Statistical Modeling` class at [Grinnell College](https://www.grinnell.edu/) during Fall 2020. 

## About SIR Model

In the basic SIR model, as susceptible members (S) of the population are exposed to a pathogen via an infected individual (I), they are moved to the infectious compartment by an infection rate, and to the recovered/removed compartment (R) after recovery/death by a recovery rate. Many simplifying assumptions are made in developing this model that could change the spread rate of an infectious disease.

The SIR model is based on differential equations, where the fraction of the population in each compartment changes as a function of time. Thus, the functions can be altered depending on the specific disease mechanism – incubation period of the pathogen, the mode of transmission, etc.

The core formulas to predict each compartment S, I, and R over time are as below:

<p align="center"><img src="https://www.lewuathe.com/assets/img/posts/2020-03-11-covid-19-dynamics-with-sir-model/ode.png"/></p>

The parameters used in these differential equations, namely the transmission and recovery rates (β & γ) are estimated based upon several factors, including but not limited to the migration rates of the population within and between countries, and the biological mechanism through which the disease affects individuals (Tang et al. 2020). For diseases like COVID-19, where the disease can be both symptomatic and asymptomatic, the models based on SIR can be extended to involve more compartments, levels to capture the complexity of transmission (Chen et al. 2020).

## Author & Contacts
[Linh Tang](https://linhtang.me/)

[Britney He](mailto:hejiayu@grinnell.edu)

[Bowen Mince](mailto:mincebow@grinnell.edu)
## Contributing
The source code for models can be found in this repository. Please note that the project is still in progress of development.

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## References
[The SIR Model for Spread of Disease - The Differential Equation Model](https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model)

[A mathematical model for simulating the phase-based transmissibility of a novel coronavirus](https://doi.org/10.1186/s40249-020-00640-3)

[A Review of Multi-Compartment Infectious Disease Models](https://doi.org/10.1111/insr.12402)

[COVID-19 Data in the World](https://www.kaggle.com/imdevskp/corona-virus-report?select=covid_19_clean_complete.csv)

[COVID-19 Data in the US](https://www.kaggle.com/imdevskp/corona-virus-report?select=usa_county_wise.csv)
