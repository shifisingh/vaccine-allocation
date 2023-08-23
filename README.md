# COVID-19 Vaccine Allocation Optimization

## Acknowledgements
This code was written by Shefali Singh and Lisa Yamada. We are grateful for the support given to us by faculty advisor Rajmohan Rajaraman, and his assistance in modifications to the code as well as high level 
strategy in formulating a reasonable model. We are also grateful for the support of Northeastern University in allowing us to continue our research as part of an undergraduate research grant. 

## Background 
The project began as this idea of formulating the distribution of Covid-19 vaccine allocation through network flow, but morphed into a more generalizable formulation of regional distribution through models which emphasized various strategies based on allocation policies and age groups as risk classes. Two papers in particular were studied in order to discern which provided a more appropriate model for abstraction. The first paper utilized the SAPHIRE model which captured varying levels of interactions across age groups, as well as implementing these interactions for both dynamic and static vaccine allocation policies. The SAPHIRE model emphasized estimation of mixing between age groups as an important aspect in discerning which allocation policy was best. The second paper utilized the DELPHI-V-OPT model. This model formulated an optimal number of vaccines to be allocated to each region and risk class to minimize the total deaths as a result of the pandemic. The DELPHI-V-OPT did not necessarily require the risk classes to be age groups. Ultimately, the SIMULATE method from DELPHI-V-OPT model was chosen as a starting off point for the model that was formulated for our research.

## Implementation 
Utilizing the standard SEIR model (susceptible, exposed, infected, recovered) and data for the state of Massachusetts, the code for our model was implemented in Python. Risk classes were considered to be the age groups of (0-19, 20-29, 30-49, 50-69, 70+). The data on vaccines allocated to these age groups came from the mass.gov weekly reports which included a weekly breakdown of vaccines allocated to age groups beginning in mid December. The numbers of initial susceptible people came from Massachusetts census data. Various rates such as government intervention (gamma) or death rate (rDeath) are estimated utilizing existing data and linear regression. The code aims to take day 0 data and with the rates estimations and vaccine allocation strategy, predict how the groups of susceptible, exposed, infected, undetected, hospitalized, quarantined, recovered, died, and immune will vary over time. Then for a time period where data is known, there is a comparison made between what the code predicts and what occurred in real life. This allows for various parameters to be adjusted and modified so that the code is able to more accurately predict for the future periods where no data is available.

## Future Work 
Currently, the code supports reading data from .csv files containing weekly data for a particular risk class and evaluating the state changes for a given time. Certain parameters still require adjustment which will come from the formal comparisons which need to be drawn between the code and the actual data. One of the project’s areas for expansion is incorporating a state “P” which captures when an individual has received one dose of a double dose vaccine. This is an area where previous research has fallen short and should hold some interesting conclusions. Another area is accounting for the looping nature of infectious diseases; a person who is considered recovered stays there for a finite time before going back to susceptible and thus restarting the loop. Though this project utilizes Covid-19 data, the hope is that a model is developed which can be applied to a variety of infectious diseases in order to mitigate the spread and death toll of the disease.

## Data References
Paper used for model 
https://www.medrxiv.org/content/10.1101/2020.11.17.20233213v1.full

Raw Historical Data
https://docs.google.com/spreadsheets/d/1J8IYQ8xkOT7k98vOWUruJmctgPAc8ZIZ/edit#gid=1127812055

Initial class values:
https://docs.google.com/presentation/d/1TbY1qXgweGtDtkhLo03vzwuD1gjDGZiLXO2bb63LfOU/edit#slide=id.gc5ebade372_0_11

Risk class info for calculating death rates/nominal infection rate/vaccination
https://docs.google.com/spreadsheets/d/1wTsipggtfNr_dDOqcH3BWIASQkVQGeLcpxEDH8S3Dg4/edit#gid=0

Vaccination Strategy (number of vax per age group per phase)
https://docs.google.com/spreadsheets/d/17VA_OWgXXKQkIFC1Tpn6Q9pWUCGQt6wl_JXyCXcZYfk/edit#gid=846397535
