# Introduction
The purpose of *dsph_search* is to use automate the search for dSph candidates using GAIA data. In order to be considered a candidate, all tests must evaluate to true. The search currently employs two tests for candidacy (in *the_search/cuts.py*):
- Angular density test: The idea here is that 'empty space' will have a roughly fixed angular density of objects/area regardless of the region considered. Gravitationally bound objects will have a higher density, so comparing the densities of a larger cone and a smaller cone should indicate if their is a gravtitationally bound (observable by GAIA) object in that area. Parameters are: density_tolerance (thresholed for ratio of the small and large cones' densities) and the radii (the two angular radii of the cones to compare).
- Proper motion test: A gravitationally bound object and its member objects should be moving with roughly the same proper motion. In the proper motion space, this should be indicated by a smeared peak. This test checks to see if the peak is 'large' enough (described below) to be considered a dSph candidate. Parameters are: radius (fixed to 0.5 but can be changed), test_area (the area on the proper motion histogram that's considered the peak; area=(2*test_area)^2), test_percentage (how 'peaked' the objects are around a central value; smaller percentage is less peaked), num_maxima (the number of maxima on the histogram tested as peaks)

# Tuning Parameters
The files in *the_search/tuning* along with *tune.py* provide a way to optimize the parameters described above. The parameters are tuned using simulated annealing and an initial guess then allowed to walk the parameter space looking for minima. The actual tuning of parameters happens in *tune.py*. The test dwarfs are loaded from *tuning/tuning_known_dwarfs.txt* and *tuning/tuning_random.txt* then the annealing procedure is run on them. The parameters are returned.

The goodness of these parameters can be checked by using the *look_at_tuned_parameter_values()* function in *main.py*. 

# Other Information
Primarily, this search is meant to be run automatically. The *main()* function in *main.py* does this conveniently. You feed it the number of samples to consider, it then divides these up as nicely as it can on a sphere and performs tests on each point. 

Successful samples are 'accepted' and they are stored with relevant plots in the *candidates/* directory. If samples are unsuccessful, they are 'rejected' and deleted; their logs are stored in *dead_logs/*.

After the full run, the successful candidates are also saved in *candidate_coords.txt* for convenience.

# Room for Improvement
- Implement more cuts
- Refine parameter values
- Automatically search Simbad for candidates and include the relevant output in the logs. 
