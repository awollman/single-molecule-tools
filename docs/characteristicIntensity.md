# Characteristic Intensity of a single fluorophore

For stoichiometry determination, the intensity of a single fluorophore must be measured. There are broadly 3 ways:

1. Image surface immobilised fluorophores in vitro, track and use the average found Intensity.
2. Totally photobleach in vivo data, track and plot the distribution of found Intensity.
3. Calculate the step-size in vivo. Run tracked or overtracked data through the Chung-Kennedy filter, calculate the pairwise distance distribution and the Fourier spectrum of this.

The most robust approach is all of them. 2 is quickest but should be accompanied by some of 3. At the very least, single photobleach steps must be observed or the data may not be single-molecule!

## Rough guide

### First Estimate
1. Generate some single molecule data and run through tracker.
2. If using multiple sets, concatenate the data using something like \Stoichiometry\getAllSpots.m
3. Plot the distribution of intensity using a KDE or histogram (column 5 of spot array). The peak value should be the characteristic intensity. Plotting this distribution of spots found only towards the very end of the bleach may yield a clearer peak.

### Overtracking
This is a method to see single photobleach steps by tracking them beyond their apparent photobleaching. Use \Stoichiometry\overTracker.m to do this and plot the resultant individual intensity traces (or wait until the next CK filtering step). Plotting the intensity distribution here again is useful as you should have 2 clear peaks at the characteristic intensity and close to 0.

### Filtering
Intensity is very noisy so mean filtering would remove any steps. So called edge-preserving filters can smooth the data without removing steps. CKall will run a Chung-Kennedy filter over the data and plot individual photobleach steps or provide input for the next section.

### Fourier spectra
Periodicity in the filtered intensity (steps!) can be found by Fourier analysis. Use `powerSpectrum` on the filterI data from CKall to calculate the pairwise distance distribution and the resultant Fourier spectrum. There should be periodic peaks in the pairwise distance distribution at the characteristic intensity and a peak at this value in the Fourier spectrum. This is often difficult to obtain as there are many sources of noise which 'dephase' the periodicity. Less is often more with this type of analysis.
