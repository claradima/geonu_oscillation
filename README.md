Hello!

The purpose of this repo is mainly to generate geoneutrino 
spectra that we expect to see given different models, mostly 
varying the Earth model and the oscillation model. The main 
script is generate_geonu_spectra.py, but a fewother scripts for
some tangential calculations will be added as well.

The **generate_geonu_spectra.py** script generates spectra for the 
following situations:

  1) Full oscillation formula, 'standard' oscillation parameters
     (i.e. theta_12 and delta_m_21^2 taken from James/Tony's
     analysis (constrained results), and the other oscillation
     parameters are fixed from PDG 2020
     
  2) Survival probability approximated with a constant (there
     will soon be a separate script that computes the average
     P_ee and error for different Earth models)

  3) Full osc formula, all oscillations parameters standard
     except theta_12 - produce two spectra, one with
     theta_12 + err from James/Tony's constrained fit, one with
     theta_12 - err

  4) Same as above, but delta_m_21^2 changes instead of theta_12

These are all saved as plots, as well as in csv files that 
contain the data displayed in the plots. In addition to the 
spectra, the script also calculates ratios between different
spectra to check how much changing various bits about the
oscillation model affects our prediction. Ratios are calculated
between any spectrum generated that is not (1) above and (1). 
These are again saved as both plots and in csv files.

There are various parameters that you can adjust when running the script, you can check what these are by opening the file and looking 
at the parser.add_argument lines. They all have default values as well. A command with custom parameters would look something like:


python3 generate_geonu_spectra.py -abd 'mid' -Ebins 100 -cgridcount 640 -mgridcount 160 -livetime 100


The other params probably don't need to be changed very often. 
You can also specify if you want the plots saved and/or shown 
while running the code.

All the plots and csv files are saved in a new subfolder in 
/plots. The naming scheme for the folder is:


dir_name = f"{abd_set}_{len(energy_array)}E{int(np.floor(grid_1d_size_crust))}C{int(np.floor(grid_1d_size_mantle))}M"


so this includes the abundance set name, the number of energy 
bins, the number of points in the crust in 1d, then the number 
of points in 1d in the mantle (these 1d_size things come from 
how the Earth model is generated : it's made as a 3d cube with 
1d_size points in 1d, then each layer is cropped as a spherical 
shell with limits on the radius)


**functions.py** contains all the functions, in a slightly 
disorganized way but they all have very detailed comments about
what they do so check that out if needed :)

Also check my explanation document for further info on this, 
or email me (c.dima@sussex.ac.uk).
