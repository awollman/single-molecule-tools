function [frame_average] = FrameAverage3(image_data, noFrames, firstLeft,firstRight,ALEX)

% Calculates a frame average of image_data array starting at startFrame
% over noFrames
%Start frame should mostly be firstLeft from LaserOn
if ALEX==0
    startFrame=min([firstLeft,firstRight]);
    frame_average=mean(image_data(:,:,startFrame:(startFrame+noFrames)),3);
    frame_average=mat2gray(frame_average);
%         frame_average=mean(image_data(:,:,firstRight:(firstRight+noFrames)),3);
%     frame_average=mat2gray(frame_average);
else
    frame_average(:,:,1)=mean(image_data(:,:,firstLeft:2:(firstLeft+2*noFrames)),3);
    frame_average(:,:,2)=mean(image_data(:,:,firstRight:2:(firstRight+2*noFrames)),3);
    frame_average=mat2gray(frame_average);
end
end