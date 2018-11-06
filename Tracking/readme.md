# Description

The tracking software consists of a series of functions for opening data and tracking particles.

# Input

The filename (or folder name for unstacked tiffs), either without the file extention or if using the bioformats plugin with extention.

# Output

SpotsCh1/2 is an array, each row contains the information for a spot found in an image frame in the series. The columns contain the following information:
1.	X coordinate (pixels)
2.	Y coordinate (pixels)
3.	Clipping_flag (a switch, please  ignore)
4.	Mean local background pixel intensity
5.	Total spot intensity, background corrected
6.	X spot sigma width
7.	Y spot sigma width
8.	Spot central intensity (peak intensity in a fitted Gaussian)
9.	Frame number the spot was found in
10.	Trajectory number, spots in the same trajectory have the same trajectory number
11.	Signal to noise ratio
12.	Frame in which laser exposure began

Useful to know: There is a ‘show_output’ option that can be used to view graphs and manually advance at each stage of the algorithm.
Cursor_mode: set =1 and user can manually specify where spots are. The code will return intensity values at this point over the whole time series.
tracker: This is the main tracking program and is a function of image_label, readData=1 to extract tif or 0 if data preloaded then runs on image_label, p is the parameter structure which can be read in or set by the code. It returns spot arrays for each channel and frame_average. It uses:
extractImageSequence3: extracts user set frames from tif specified by image_label, can open Andor ASCII files and folders full of tif frames
ImEx1: uses bftools to open many life science image formats. See Open Microscopy Environment for details.
LaserOn3: Calculates where the first illuminated frame is based on maximum intensity
FrameAverage2: Calculates a frame average/summation over set no. frames
findSpots2: thresholds the image with Otsu’s method to find candidate spots
findSpots3: thresholds the image with Otsu’s method to find candidate spots
findSpots4: performs findSpots2+3 in a more efficient way
findSpotCentre3: performs iterative Gaussian masking to find spot centre and total intensity
fit2DgaussianFixedCenter2: fits a constrained 2D Gaussian to find sigma_x/y and central intensity
MergeCoincidentCandidates2: will use pairwise distances to remove candidates which are too close to be resolved, in some circumstances faster by quite a slow function in itself
iterate1DgaussianFixedCenter2: finds PSF width by masking with Gaussians of different sizes
LinkSpots4: links spots into trajectories based on proximity
