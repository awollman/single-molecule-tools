% Function to apply a non-linear transformation between 2 channels of a
% microscope based on manual point selection

%imageCh1/2 are single frame images from each channel to align, preferably
%of a bead sample which bleeds into both channels
%SpotsCh1/2 are the track arrays for each channel ALREADY CORRECTED FOR THE
%OFFSET BETWEEN IMAGECH1/2. i.e. imageSize subtracted from SpotsCh2.


function alignedSpots=spotChManualReg(imageCh1, imageCh2, SpotsCh1, SpotsCh2, showOutput)
[moving_out,fixed_out] = cpselect(imageCh2,imageCh1,'Wait', true);
tform = fitgeotrans(moving_out,fixed_out,'NonreflectiveSimilarity');
Jregistered = imwarp(imageCh2,tform,'OutputView',imref2d(size(imageCh1)));
falsecolorOverlay = imfuse(imageCh1,Jregistered);
[SpotsCh2(:,1), SpotsCh2(:,2)] = transformPointsForward(tform, SpotsCh2(:,1), SpotsCh2(:,2));
alignedSpots=SpotsCh2;
if showOutput==1
    figure;
    imshow(falsecolorOverlay,[])
    figure;
    subplot(1,2,1)
    scatter(SpotsCh1(:,1),SpotsCh1(:,2),'g')
    hold on
    scatter(SpotsCh2(:,1),SpotsCh2(:,2),'r')
    title('misaligned spots')
    
    subplot(1,2,2)
    scatter(SpotsCh1(:,1),SpotsCh1(:,2),'g')
    hold on
    scatter(alignedSpots(:,1),alignedSpots(:,2),'r')
    title('aligned spots')
    
end

end