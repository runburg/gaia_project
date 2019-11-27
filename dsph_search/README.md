# Introduction
The purpose of *dsph_search* is to use automate the search for dSph candidates using GAIA data. In order to be considered a candidate, all tests must evaluate to true. The search currently employs two tests for candidacy (in *the_search/cuts.py*):
- Angular density test: The idea here is that 'empty space' will have a roughly fixed angular density of objects/area regardless of the region considered. Gravitationally bound objects will have a higher density, so comparing the densities of a larger cone and a smaller cone should indicate if their is a gravtitationally bound (observable by GAIA) object in that area. Parameters are: density_tolerance (threshold for ratio of the small and large cones' densities) and the radii (the two angular radii of the cones to compare).
- Proper motion test: A gravitationally bound object and its member objects should be moving with roughly the same proper motion. In the proper motion space, this should be indicated by a smeared peak. This test checks to see if the peak is 'large' enough (described below) to be considered a dSph candidate. Parameters are: radius (fixed to 0.5 but can be changed), test_area (the area on the proper motion histogram that's considered the peak; area=(2&ast;test_area)^2), test_percentage (how 'peaked' the objects are around a central value; smaller percentage is less peaked), num_maxima (the number of maxima on the histogram tested as peaks)

# Tuning Parameters
The files in *the_search/tuning* along with *tune.py* provide a way to optimize the parameters described above. The parameters are tuned using simulated annealing and an initial guess then allowed to walk the parameter space looking for minima. The actual tuning of parameters happens in *tune.py*. The test dwarfs are loaded from *tuning/tuning_known_dwarfs.txt* and *tuning/tuning_random.txt* then the annealing procedure is run on them. The parameters are returned.

The goodness of these parameters can be checked by using the *look_at_tuned_parameter_values()* function in *main.py*. 

# Other Information
Primarily, this search is meant to be run automatically. The *main()* function in *main.py* does this conveniently. You feed it the number of samples to consider, it then divides these up as nicely as it can on a sphere and performs tests on each point. 

Successful samples are 'accepted' and they are stored with relevant plots in the *candidates/* directory. If samples are unsuccessful, they are 'rejected' and deleted; their logs are stored in *dead_logs/*.

After the full run, the successful candidates are also saved in *candidate_coords.txt* for convenience.

# A More Detailed Account
The algorithm consists of two parts: querying GAIA and applying the tests.

1. **Querying GAIA**: This step consists of getting GAIA data from the GAIA DR2 archive using the astroquery asynchronous job TAP query. The script for this search is found in *the_search/utils.py* and is called with the *gaia_search()* function.     The query is written in ADQL and is given as:
    ```
    SELECT TOP 500000
    gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,
    gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,
    gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,
    gaia_source.bp_rp
    FROM gaiadr2.gaia_source
    WHERE
    CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',{coords.ra.degree}, {coords.dec.degree}, {radius}))=1 
    AND (gaiadr2.gaia_source.parallax - gaiadr2.gaia_source.parallax_error * {sigma} <= 0) 
    AND (SQRT(POWER(gaiadr2.gaia_source.pmra, 2) + POWER(gaiadr2.gaia_source.pmdec, 2)) <= {pm_threshold}) 
    AND (gaiadr2.gaia_source.bp_rp <= {bp_rp_threshold})
    ```
    The clauses after **WHERE** determine which objects are returned by the query. The objects must be 
    - Within a circle of {radius} centered at ({coords.ra.degree), {coords.ra.degree})
    - Have parallax consistent with zero (default value of {sigma} is 5)
    - Have proper motion magnitude less than {pm_threshold} (default value 5 mas/yr)
    - Have bp_rp less than {bp_rp_threshold} (default value is 2)
    Note: the {} are Python syntax, not part of the ADQL.
    
    The returned VOTable table from the query is stored in the dwarf's directory in the *vots/* subdirectory. This step is taken because the GAIA query is the most time-consuming part of the algorithm. The tables are saved to expedite later computational steps.
  
2. **Running the tests**: As stated above, there are currently two tests run on each set of coordinates. These tests will first attempt to load a previously stored GAIA query from *vots/* before initiating a new GAIA query to save time.

      - **Angular density test**: The angular density test compares the angular density (number or objects/radius of cone^2) for two different cones centered on the coordinates by comparing the ratio of these densities to a threshold value. We expect patches of sky with no gravitationally bound objects to have constant angular density and thus a ratio ~1. Parameters are: density_tolerance (threshold for ratio of the small and large cones' densities) and the radii (the two angular radii of the cones to compare; default values are 1.5, 0.1). 
      - **Proper motion test**: The proper motion test evaluates the morphology of a 2d proper motion histogram. The test identifies the maxima of the the histogram (indicating many objects with similar proper motion) and then looks at the area on the histogram immediately surrounding the maxima to see if this maxima is a 'peak' in the proper motion space. 'Peak' means that there is a sufficiently large percentage of the total objects on the histogram near this maxima value. Parameters are: radius (fixed to 0.5 but can be changed), test_area (the size of the area around the maxima considered; area=(2&ast;test_area)^2), test_percentage (how 'peaked' the objects are around a central value; smaller percentage is less peaked), num_maxima (the number of maxima on the histogram tested as peaks)

    
# Room for Improvement
- Implement more cuts
- Refine parameter values
- Automatically search Simbad for candidates and include the relevant output in the logs. 

# Change in Methodology
This branch changes and optimizes the approach used to search for dSphs. Now, the process looks like:
 1. GAIA data for a region is downloaded.
 2. Within this region, the GAIA objects are binned for the entire region using a sufficiently small bin size (0.02 deg usually).
 3. The histogram is convolved with a tophat of varying radial size (can also convolve with a 2d Gaussian kernel, but it is more computationally expensive).
 4. The convolved histogram is searched for overdensities (based on Poisson statistics, num_sigma=2). If the same bin has an overdensity in a few (repitition=2) of the convolved histograms, it is considered a dwarf candidate. Both of these parameters can be adjusted.
    
**Some initial observations:**
- Nearby bright objects frequently prevent detection, so faint objects favor a smaller test region radius to be detected (HyiI).
- The stripiness of GAIA might mess with statistics (LeoT).
- Many of the known dwarfs seem to be near a bright object (RetII, Scl).
- I'm not sure what's happening with SgrII.
