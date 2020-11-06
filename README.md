# enzyme_screen


Scripts and functions for extracting and analysing biochemical reactions.

Author: Andrew Tarzia
Email: andrew.tarzia@gmail.com or atarzia@ic.ac.uk

This work was produced in the final year of my PhD at the University of Adelaide under the supervision of A/Prof David Huang and Prof Christian Doonan.

Previously at: https://bitbucket.org/andrewtarzia/psp_source/src/master/

A Jupyter notebook that runs through the molecular size calculation from a SMILES string is available:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/andrewtarzia/enzyme_screen/master?filepath=examples%2Fcalculate_molecular_size.ipynb)

## Installation

* Tested on Ubuntu 18.04 using conda and pip
* Install Anaconda in standard way (Python 3.7.3)
* Packages required outside of what comes with conda
 * RDKit:
   * `conda install -c conda-forge rdkit`
   * Version: 2019.09.2.0
 * chemcost:
   * Python code written by Steven Bennett for the extraction of purchasability from the ZINC15 database.
   * Follow instructions found here: https://github.com/stevenbennett96/chemcost
   * Only required for `molecule_population.py`

## Workflow


### Collecting database from KEGG

* Download br08201 JSON file from [the KEGG library](https://www.genome.jp/kegg-bin/get_htext?query=08201&htext=br08902.keg)
 * Used version as of May12_2020 of br08201: Enyzmatic reactions
* Run `util/split_KEGG.py` in working directory to produce:
    * `_ECtop.json`: A dictionary of all reactions for all ECs
    * `_EClist.txt`: A list of all ECs to iterating through
* Update `data/param_file.txt` with location of these files.


### Parameter testing


* All parameter screens in the supporting information of DOI: **awaiting** are run in `param_screening.py`
    * `data/test_molecules.txt` contains the required molecular information
    * Within `param_screening.py` are the range of parameters to test, the originals being set in `data/param_file.txt`


### Reaction collection and analysis

* `RS_collection.py`
    * Iterates through provided EC and reaction files to collect reaction systems
    * Also collects unique molecules to molecule database
    * Currently only implements API for KEGG
    * To be run in directory with reactions

* `molecule_population.py`
    * Trivial parallelisation done using `utils/molecule_splitter.py`
    * Takes _unopt.mol file of all collected molecules:
        * Optimises them using ETKDG -> _opt.mol
        * Calculates their properties -> _prop.json
        * Calculates the molecule size of N conformers -> _size.csv
    * To be run in directory with molecules
    * Produces some plots of chemical space

* `chemical_space_plot.py`
    * Iterates through all collected molecules and plots various chemical space plots
    * To be run in directory with molecules

* `RS_analysis.py`
    * `molecule_population.py` must be run before this point!
        * Unanalysed molecules result in skipped reactions
    * Populates the properties of each reaction system based on the properties of constituent molecules (in molecule database)
    * To be run in directory with reactions
    * Outputs all properties to `rs_properties.csv`

* `screening.py`
    * Produces the plots and screening of all reaction systems seen in DOI: **awaiting**
    * Multiple cases are defined within the script to look at specific EC numbers or system types
        * case = production for plots in DOI:
    * To be run in directory with reactions


## Examples

* `biomin_screening.py`
    * A script used to produce Figure **XX** in DOI:
    * Analyses a list of molecules that have been tested for enzyme@ZIF-8 reactions

* `examples/calculate_molecular_size.ipynb`
    * Jupyter notebook that runs a user through calculating the size of any molecule

* `examples/screen_new_reactions.ipynb`
    * Jupyter notebook that runs through the screening process exemplified in the paper search for new reactions

* `visualise_ellipsoid_steps.py`
    * Allows the user to visualise the step-wise calculation of the min. vol. enclosing ellipsoid

* `visualise_reaction_system.py`
    * Allows the user to print properties of a reaction system
