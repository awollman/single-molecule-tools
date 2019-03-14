% Function to apply a non-linear transformation between 2 channels of a
% microscope based on similarity registration of bead or brightfield images

% ALH 26.02.2019

%imageCh1/2 are single frame images from each channel to align, preferably
%of a bead sample which bleeds into both channels
%SpotsCh1/2 are the track arrays for each channel ALREADY CORRECTED FOR THE
%OFFSET BETWEEN IMAGECH1/2. i.e. imageSize subtracted from SpotsCh2.


function [alignedSpots,tform]=spotChSimilarity(imageCh1, imageCh2, spots1, spots2, showOutput)

[optimizer, metric] = imregconfig('multimodal');

optimizer = registration.optimizer.OnePlusOneEvolutionary;
optimizer.InitialRadius = 0.001;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;

metric = registration.metric.MattesMutualInformation;
metric.NumberOfSpatialSamples = 500;
metric.NumberOfHistogramBins = 50;
metric.UseAllPixels = 1;

transformtype ='similarity';  %'translation';'rigid';'similarity';'affine';

tform = imregtform(imageCh2,imageCh1,transformtype,optimizer,metric);    
Jregistered = imwarp(imageCh2,tform,'OutputView',imref2d(size(imageCh1)));
falsecolorOverlay = imfuse(imageCh1,Jregistered);
alignedSpots=spots2;
[alignedSpots(:,1), alignedSpots(:,2)] = transformPointsForward(tform, spots2(:,1), spots2(:,2));

if showOutput==1
    figure;
    imshow(falsecolorOverlay,[])
    figure;
    subplot(1,2,1)
    scatter(spots1(:,1),spots1(:,2),'g')
    hold on
    scatter(spots2(:,1),spots2(:,2),'r')
    title('misaligned spots')
    
    subplot(1,2,2)
    scatter(spots1(:,1),spots1(:,2),'g')
    hold on
    scatter(alignedSpots(:,1),alignedSpots(:,2),'r')
    title('aligned spots')
    
end

end

%WORKED EXAMPLE: TO USE you must e.g. load in brightfield image, split in half to make imagech1, imagech2.  
%image_data = rot90(image_data,3);  %rotation may be necessary for horizontal camera slits
%[image_Y,image_X,numFrames]=size(image_data);
%framewidth=round(image_X/2);
%framech1=image_data(:,1:framewidth,1);
%framech2=image_data(:,framewidth+1:end,1);
%figure; subplot(1,3,1); imshow(framech1,[]); subplot(1,3,2); imshow(framech2,[]);subplot(1,3,3);imshowpair(framech1,framech2);

%Then load in TRACKS.mat file and move channel 2 spots to overlap the channel 1 spots:
%SpotsCh2(:,1)=SpotsCh2(:,1)-framewidth;
%figure; scatter(SpotsCh1(:,1),SpotsCh1(:,2)); hold on; scatter(SpotsCh2(:,1)+0,SpotsCh2(:,2)+0);

%Finally perform the image registration and correct the channel 2 spot positions:
%[aligned_spotsCh2,tform] =spotChSimilarity(framech1,framech2,SpotsCh1,SpotsCh2,1);
%SpotsCh2 = aligned_spotsCh2;
%SpotsCh2(:,1)=SpotsCh2(:,1)+framewidth;