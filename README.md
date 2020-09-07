# FreqMagSEA v0.1
Python code to derive frequency-magnitude relationships for volcanoes in Southeast Asia using Bayesian inference and model averaging. Code accompanies a forthcoming paper currently in preparation.

## Table of contents
* [General info](#general-info)
* [Libraries/packages](#libraries-packages)
* [Set-up and Run](#Set-up)
* [Example](#example)
* [Features](#features)
* [Status](#status)
* [Contact](#contact)

## General info
Add more general information about project. What the purpose of the project is? Motivation?

## Libraries/packages
* [Pandas](https://pandas.pydata.org/) - Version 0.25.1
* [Numpy](https://numpy.org/) - Version 1.16.5
* [MatPlotLib](https://matplotlib.org/) - Version 3.1.1
* [PyMC3](https://docs.pymc.io/) - Version 3.8
* [Seaborn](https://seaborn.pydata.org/) - Version 0.9.0
* [SciPy](https://www.scipy.org/) - Version 1.3.1

## Set-up and Run

### Data
A number of CSV files are provided with this code to conduct the analysis. Details of these are:  

* Volcano_list.csv - Downloaded from Global Volcanism Program (GVP) on 04 Mar 2020
* GVP\_DB_2019 - Downloaded from GVP on 27 Jan 2020 
* Change_points.csv - Change points (50th percentile) from [Mead and Magill 2014](https://doi.org/10.1007/s00445-014-0874-y).
* Whelley_SEA - Volcanoes and classification of volcano type from [Whelley et al. 2015](https://doi.org/10.1007/s00445-014-0893-8)

#### If you want to use updated volcano records
Download the [Holocene volcanoes list](https://volcano.si.edu/list_volcano_holocene.cfm) and [eruptions database](https://volcano.si.edu/search_eruption.cfm) from GVP website. At present, the analysis will only run for volcanoes that are included within both the Whelley and GVP volcano list (n=173). 
  
**Note:** Some modification (e.g. headers) to spreadsheets downloaded from GVP may be necessary. See example spreadsheets: "GVP\_DB\_2019.csv" and "Volcano_list.csv" for formatting guidance. Ensure that spreadsheets are named consistently in two python scripts. 

### To run a single volcano
1. Copy the volcano's six digit Volcano Number from [GVP](https://volcano.si.edu/) into the "Volc_Num.txt" file.
2. Run "FreqMagSEA.py". 

If a Bayesian update is conducted the code will take approximately 5-10 minutes to run. Otherwise it will take a matter of seconds as it will only assign a prior for volcanoes with no eruption record in the GVP.

### To run all applicable Southeast Asian volcanoes
1. Run "RunMultiple.py"

**Note:** It will take a long time to run through the entire set of volcanoes (approx 10-11 hours). This may change (increase) in future if more volcanoes are able to have the Bayesian update conducted (e.g. eruptions from geological record added to GVP database or new eruptions at previously dormant volcanoes with no recorded eruptions). It is possible the time could be reduced by changing the number of cores used in the ""FreqMagSEA.py", but [PyMC3/Theano has known issues](https://github.com/pymc-devs/pymc3/issues/3140#issuecomment-453850707) with multi-processing on Windows, so this is not currently recommended. 

Volcanoes that produce an error will be printed in the "Volcanoes\_not_run.txt" file. With the existing set of volcanoes, this will be because they're classified in the GVP as "submarine", which we do not use in this analysis. 

## Example
We provide a worked example using Merapi volcano as a case study that presents the general scientific rationale and assumptions of the approach. This worked example can be found in the associated Jupyter notebook.

## Features
* Automatically assigns a prior to volcanoes based on analogues classifications.  
* Conduct a Bayesian update for volcanoes with a volcanic record.  
* Model averaging for two different approaches to analogue volcano classification.

To do list:  

* Automatic change point detection

## Status
Project is: in progress

## Contact
Josh Hayes  
Susanna Jenkins