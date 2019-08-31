# Gaia project
# Search routine for dSphs

22-08-2019 13:58
The goal of this code will be to search for and identify new dwarf spheroidal galaxy candidates using GAIA data.
This log will serve to track tasks and progress along the way.

Project outline:
  - *main.py* will be used for implementing the primary testing routine
  - dSph candidates with relevant parameters and plots will be stored in *candidates/*
  - scripts and other operational files will be in *the_search/*


24-08-2019 15:08
  - each dwarf now generates its own directory in *candidates/*
    - there information re: dwarf including tables and plots can be stored
  - *Dwarf.rejected()* deletes this directory and keeps the log in *candidates/dead_logs/*
  - *utils.py* provides accessory functions
  - *cuts.py* provides the scripts for the actual cut tests
  - *plotting_functions.py* provide plotting functions for successful dwarfs


28-08-2019 15:42
 - sample dwarf files have been generated for a set of random cones and known dwarfs
 - *cuts.py* provides the following test:
  - *parallax_test()* tests for dwarf candidacy by only looking at objects consistent (to different degrees) with zero
    - N.B.: It is not working.
  - *proper_motion_test()* tests for dwarf candidacy by analyzing the object in proper motion parameter space


  29-08-2019 14:23
