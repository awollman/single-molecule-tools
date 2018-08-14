# Copy number

Follow this guide to deconvolute single-molecule images into actual copy numbers of fluorescent molecules.

This works by generating a simulated image (or images) of the fluorescent object by integrating a 3D point spread function (PSF) over a 3D model of the object. Then solving a set of linear equations to calculate the actual number of fluorophores generating that object.

## Point Spread function

Generate a measured or simulated PSF. This [code](https://uk.mathworks.com/matlabcentral/fileexchange/31945-widefield-fluorescence-microscope-point-spread-function) will generate a good PSF. Make sure each Z-plane is normalised to unitary intensity and you simulate enough z-planes for your 3D object.

## Simulate shape
[This folder](../Image Simulation/3Dshapes) contains many functions for simulating shapes like spheres, shells and cylinders and combining them to create more complicated shapes like rod shaped bacterial cells.

## Single compartment copy number
