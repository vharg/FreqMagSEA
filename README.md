# FreqMagSEA v1
Python code to derive frequency-magnitude relationships for volcanoes in Southeast Asia using Bayesian inference and model averaging. Code accompanies a forthcoming paper currently in preparation.

## Table of contents
* [General info](#general-info)
* [Libraries/packages](#libraries-packages)
* [Set-up and Run](#Set-up)
* [Example](#example)
* [Features](#features)
* [Status](#status)
* [Contact](#contact)

## Libraries/packages
* [Pandas](https://pandas.pydata.org/) - Version 1.2.3
* [Numpy](https://numpy.org/) - Version 1.20.1
* [MatPlotLib](https://matplotlib.org/) - Version 3.3.4
* [PyMC3](https://docs.pymc.io/) - Version 3.11.1
* [Seaborn](https://seaborn.pydata.org/) - Version 0.11.1
* [SciPy](https://www.scipy.org/) - Version 1.6.0

## Set-up and Run

### Data
A number of CSV files are provided with this code to conduct the analysis. Details of these are:  

* Volcano_list.csv - Downloaded from Global Volcanism Program (GVP) on 04 Mar 2020
* GVP\_DB_2019 - Downloaded from GVP on 27 Jan 2020
* Change_points_all.csv - Change points from [Mead and Magill 2014](https://doi.org/10.1007/s00445-014-0874-y).
* Whelley_SEA - Volcanoes and classification of volcano type from [Whelley et al. 2015](https://doi.org/10.1007/s00445-014-0893-8)

#### If you want to use updated volcano records
Download the [Holocene volcanoes list](https://volcano.si.edu/list_volcano_holocene.cfm) and [eruptions database](https://volcano.si.edu/search_eruption.cfm) from GVP website. At present, the analysis will only run for volcanoes that are included within both the Whelley and GVP volcano list (n=173).

**Note:** Some modification (e.g. headers) to spreadsheets downloaded from GVP may be necessary. See example spreadsheets: "GVP\_DB\_2019.csv" and "Volcano_list.csv" for formatting guidance. Ensure that spreadsheets are named consistently in two python scripts.

### To run a single volcano
1. Comment out the "selected_volcanoes" variable (line 40 of "Analysis.py")
2. Change the "volcanoes" variable (line 41 of "Analysis.py") to a list that contains just the GVP of the volcano of interest.  

If a Bayesian update is conducted the code will take approximately 5-10 minutes to run per loop. Otherwise it will take a matter of seconds as it will only assign a prior for volcanoes with no eruption record in the GVP. So if all loops are included (i.e. full sensitivity analysis) the code will take a hours per volcano.

### To run all applicable Southeast Asian volcanoes
1. Run "Analysis.py"

**Note:** It will take a long time to run through the entire set of volcanoes for all potential loops.  It is possible the time could be reduced by changing the number of cores used in the ""FreqMagSEA.py", but [PyMC3/Theano has known issues](https://github.com/pymc-devs/pymc3/issues/3140#issuecomment-453850707) with multi-processing on Windows, so this is not currently recommended.

## Example
We provide a worked example using Merapi volcano as a case study that presents the general scientific rationale and assumptions of the approach. This worked example can be found in the associated Jupyter notebook.

## Features
* Automatically assigns a prior to volcanoes based on analogues classifications.  
* Conduct a Bayesian update for volcanoes with a volcanic record.  
* Model averaging for two different approaches to analogue volcano classification.

## Status
Project is: in progress

## Contact
Josh Hayes  
Susanna Jenkins
