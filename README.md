# Watters-et-al-2020-jNeurosci
Matlab code to semi-automatically identify stationary mitochondria from time-lapse imaging of axonal mitochondria (requires user input)

If you use this code, please cite: Watters et al., (2020) AMPK preferentially depresses retrograde transport of axonal mitochondria during localised nutrient deprivation. J Neurosci XX:XX

Sample analysis:
1. Download and store files in folder
2. In MATLAB, run MitoCount('sample_cellprofile_output1') -> this generates 'sample_cellprofiler_output1.mat'
3. In MATLAB, run MitoCount_stationary('sample_cellprofile_output1') -> this draws kymograph and goes through potentially stationary objects with user. User has to confirm or reject suggestions. Final kymograph (stationary objects drawn red) is saved and the numbers of stationary mitochondria in each time interval analysed are output to the screen.
