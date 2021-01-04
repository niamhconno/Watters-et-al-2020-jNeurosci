# Watters-et-al-2020-jNeurosci 
If you use this code or documents herein, please cite: Watters et al., (2020) AMPK preferentially depresses retrograde transport of axonal mitochondria during localised nutrient deprivation. J Neurosci PMID: 32393534 DOI: 10.1523/JNEUROSCI.2067-19.2020


Matlab code to semi-automatically identify stationary mitochondria from time-lapse imaging of axonal mitochondria (requires user input)

A sample input file is provided ('sample_cellprofiler_output.csv'), and sample output files (.mat and .fig) begin with 'X...' (one .bmp file provided to open outside Matlab. Output filenames have been edited to avoid interference when you run the sample analysis on your own machine :-).

Sample analysis:
1. Download and store files in folder
2. In MATLAB, run MitoCount_stationary('sample_cellprofiler_output1',450,450,1,2) 
-> this generates 'sample_cellprofiler_output1.mat', then draws kymograph and goes through potentially stationary objects with user. User has to confirm or reject suggestions. Final kymograph (stationary objects drawn red) is saved and the numbers of stationary mitochondria in each time interval are output to the screen.

See Word document ('detailed-protocol-MitoCount-stationary.doc') for more details on the algorithm and how to run the analyses. 
