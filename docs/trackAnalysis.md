# Track Analysis

Turn tracking data into stoichiometry and diffusion coefficients. Look for colocalisation.

`trackAnalyser` - runs on a single cell containing one colour channel
Inputs: spot array, segmentation mask for one cell, image file name, cell number and a parameter structure
Outputs: trackArray containing stoichiometry, diffusion coefficient etc. and the corresponding spots

`colocalisedTrackAnalyser` - runs on a single cell containing 2 colour channels
Inputs: spot arrays, segmentation mask for one cell, image file name, cell number and a parameter structure
Outputs: trackArrays containing stoichiometry, diffusion coefficient etc. and the corresponding spots

`sampleTrackAnalyser` and `sample2CTrackAnalyser` show you how to loop the previous functions over whole data sets and then plot things with `masterPlot`
