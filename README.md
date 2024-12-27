# AF pipeline

**"The road ahead may be cloaked in shadow, but light enough there is for the one bold enough to seek it."**  

The following repository contains scripts to parse and analyze AF2/3 predictions.  

## Overview
The organization and roles of the scripts are as follows:
1. `parser.py`: defines a class AfParser that comprises methods to parse and extract relevant features (coordinates, pLDDT, and PAE for now) from the AF2/3 output files.  
Any new method required to extract a feature from the input files must go here.  

2. `main.py`: contains classes to be used for obtaining required metadata (contact maps, interacting residues) from the extracted features.  
Creates the required attributes (features) required for getting the metadata.  
Any class for extracting new metadata must be created here.  

3. `utils.py`: contains accessory functions.  

4. `gata.py`: for each project one can write a script to perform any downstream analysis given the metadata.  

Look at the scripts for more details.


## AF3 

Take 100 models from 20 seeds (5 models each). 

## AF2

Need to look at `ranking_debug.json` to get the best model. 
