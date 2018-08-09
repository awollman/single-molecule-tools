# Deconvolution

Follow this guide to deconvolute single-molecule images into actual copy numbers of fluorescent molecules.

This works by generating a simulated image (or images) of the fluorescent object by integrating a 3D point spread function (PSF) over a 3D model of the object. Then solving a set of linear equations to calculate the actual number of fluorophores generating that object.

## Point Spread function

Generate a measured or simulated PSF. This [code](https://uk.mathworks.com/matlabcentral/fileexchange/31945-widefield-fluorescence-microscope-point-spread-function) will generate a good PSF. Make sure each Z-plane is normalised to unitary intensity and the PSF array is the same size or larger than the image in all 3D.
