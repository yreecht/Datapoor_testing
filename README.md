R codes for simulation testing variety of data poor/rare species estimation method
< in preparation >

================

*Contributors: Kotaro Ono<sup>1\*</sup>*

<sup>1</sup> Institute of Marine Research, Bergen, Norway.<br>
<sup>\*</sup> Corresponding author: <kotaro.ono@hi.no><br>


The purpose of this github repository is create a spatially explicit
simulation framework to generate more "realistic" data from the fishery.
The model simulates the abundance in space of k number of species with 
each species having its own population and movement dynamics. 
Fishery activity is based on ideal free distribution where fishers
distribute their effort based on the expected revenue of the fishing grounds.
The spatial extent of each fisher can be specify so that some might only fish 
in certain area of the simulated space. 
Primary sampling unit is then either a vessel or a haul and we can specify 
different level of sampling intensity depending on the scenario examined.
In this repository, you will find the necessary R- and C++ code to run the 
simulation, and fit a variey of data poor models.

## Code

The code is available in the *R/* folder and divided into different
scripts. In short, the *R/OM.R* file specifies the scenario and runs 
estimation methods. The remaining scripts are utlity functions to 
run the spatial explicit population dynamics model under the conditions
specified in *R/OM.R* or run the probabilistic sampling for data 
collection, and more.

Some models (i.e. the multispecies DFA ) are implemented using 
*Template Model Builder* (TMB) and the C++ code can be found in the 
*src/* folder.

## Documentation

Documentation about the model will be added into the *docs/* folder

## Author’s github accounts

Kotaro Ono - [kotkot](https://github.com/kotkot)

## Licence

This project is licensed under the GNU GPLv3 License - see
[LICENSE](LICENSE) for details.

## Funding

...

## References
  
[<img src="https://www.hi.no/en/hi/resources/layout/HI-logo-farger-engelsk.svg/original"
alt="Institute of Marine Research" width="200"/>](https://www.hi.no/en)
