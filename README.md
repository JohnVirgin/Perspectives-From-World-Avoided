# Perspectives-From-World-Avoided
Collection of Python Scripts and iPython Notebooks (to be used via jupyter notebook IDE) used in the formulation of results for the Study:

Is Arctic Amplification Dominated by Regional Forcing and Feedbacks: Perspectives from the 'World-Avoided' Scenario.

# Full Citation

Virgin, J. G., & Smith, K. L. (2019). Is Arctic Amplification dominated by regional radiative forcing and feedbacks: Perspectives from the World-Avoided scenario. Geophysical Research Letters, 46(11).


# Preamble notes
- As long as all of these jupyter notebooks are in the same working directory, you shouldn't have to modify the directory paths, and the scripts should work straight away along with the data that is included.
- The .py scripts in /Variable Get Scripts folder require the entire Model output to use. See the acknowledgments section of the study's full text for more information

# Plotting Scripts
- These plotting scripts should work out-of-the-box should you also have the necessary companion folder of pickled data files (available from the Arctic Data Centre under 'Perspectives_From_WAVD_Pickled'), the necessary supplementary packages installed (see the imported packages for each script), and are at running least a python3 build.
- Note, the Plot_Fig_2 script will NOT work with exclusively the folder of pickled arrays (see above point), as it calls directly upon the radiative forcing data created with PORT in this study.

# Feedback and Misc Calculations
- Here, the scripts that are associated feedbacks derived from the GFDL and ECHAM6 radiative Kernels are NOT included, as they are near identical to the CAM3 scripts provided
- The scripts call upon pickled (https://docs.python.org/3/library/pickle.html) Climate change responses that were calculated and saved directly from the model output.

# Statistical Significance Scripts
- Same as above, only statistical significance scripts for feedbacks derived from the CAM3 kernels are included here
- Scripts call upon pickled radiative feedback lists
