function [Istep_axis, power_spectrum]=powerSpectrum(filteredI,stepSize, showOutput)
PwD=pdist(filteredI);
counts=histcounts(PwD,0:stepSize:max(PwD));
x=0:stepSize:max(PwD);
binCentres=x;
freqCounts=counts;
L = length(freqCounts); 
nFFT = 2^nextpow2(L);
spectrum_pwd_distrib = fft(freqCounts,nFFT)/L; 
bin_separation = (binCentres(end)-binCentres(1))/(length(binCentres)-1);
Fs = 1/bin_separation; 
freq_axis = Fs/2*linspace(0,1,nFFT/2+1); 
warning('off','MATLAB:divideByZero'); 
Istep_axis = 1./freq_axis; 
power_spectrum = (2*abs(spectrum_pwd_distrib(1:nFFT/2+1))).^2;
if showOutput==1
 figure
plot(Istep_axis,power_spectrum,'b-','LineWidth',1.5);
xlabel('Step Size')
ylabel('Power')
end
end