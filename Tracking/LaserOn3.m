% Created by Adam Wollman 2014
% If you use this function please include me as second author in your
% publications and send a cheque to PO box...
%
%Determines when the laser switched on in an image array (y,x,time) and
% outputs average pixel values. ALEX =1 if alternating laser experiment and
% =0 if CW.
% INPUTS
%
% image_data an array containing the image data
%
% ALEX =1 if alternating laser experiment and =0 if CW.
%
% use_diff=1 to use the max difference between frames as start frame, good
% for low signal
%
% OUTPUTS
%
% firstLeft is the first frame in the left channel
%
% firstRight is the first frame in the right channel
%
% LeftAverage is the average pixel intensity in the left channel for each frame
%
% RightAverage is the average pixel intensity in the right channel for each frame
%
%

function [firstLeft, firstRight, LeftAverage, RightAverage] = LaserOn3(image_data, use_diff, ALEX)
SizeData=size(image_data);
NoFrames=size(image_data,3);
SizeFrame=size(image_data);
LeftAverage=zeros(1,NoFrames);
RightAverage=zeros(1,NoFrames);
AllAverage=zeros(1,NoFrames);
for p=1:NoFrames
    %Average of 473 (LHS) channel
    LeftAverage(1,p)=mean(mean(image_data(:,1:SizeFrame(2)/2,p)));
    %Average of 561 (RHS) channel
    RightAverage(1,p)=mean(mean(image_data(:,SizeFrame(2)/2:SizeFrame(2),p)));
    AllAverage(1,p)=mean(mean(image_data(:,:,p)));
end

if ALEX==0
    if use_diff==0
        [~,firstLeft]=max(LeftAverage);
        [~,firstRight]=max(RightAverage);
      %  firstRight=firstLeft;
    else
        [~,firstLeft]=max(diff(AllAverage));
        if firstLeft>2
            firstLeft=firstLeft;
        end
        firstRight=firstLeft;
    end
else
    if use_diff==0
        %Assume max average intensity frame on LHS is the first
        %473 frame
        [~,firstLeft]=max(LeftAverage);
        if firstLeft>2
            firstLeft=firstLeft;
        end
        %Identify first 561 frame by checking if frame
        %before first 473 is larger than after
        if RightAverage(firstLeft-1) > LeftAverage(firstLeft+1)
            firstRight=firstLeft-1;
        else
            firstRight=firstLeft+1;
        end
    else
        [~,firstLeft]=max(diff(AllAverage));
        if firstLeft>2
            firstLeft=firstLeft;
        end
        %Identify first 561 frame by checking if frame
        %before first 473 is larger than after
        if RightAverage(firstLeft-1) > LeftAverage(firstLeft+1)
            firstRight=firstLeft-1;
        else
            firstRight=firstLeft+1;
        end
        
        
    end
    
    
end