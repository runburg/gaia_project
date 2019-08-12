# Gaia project

The goal of this project is to identify new dwarf spheroidal galaxy candidates through the use of GAIA data.
We aim to get an initial O(0) idea of where these objects could be.

The first steps will be to understand the average stellar velocity and stellar object density in known (classical) dSphs.
Then in this parameter space, use GAIA data to indicate or exclude stellar object structures as dSphs.

If the results are encouraging, the next steps would be to employ machine learning and simulation to more carefully determine dSph candidates.


The workflow is:
  - *main.py* queries SIMBAD and GAIA and gets all necessary astrophysical quantities. If desired, this can be done for multiple cone radii at once
  - *gaia_functions.py* have all of the necessary functions for the data processing and plotting
  - *gaia_plotting.py* simply calls the functions from above and deals with the file names, etc.
