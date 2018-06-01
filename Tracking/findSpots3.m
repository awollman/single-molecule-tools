function [result] = findSpots3(frame,method,disk_radius,gaussian, output)
%
% Find candidate (fluorescence) bright spots on images (of cells).
% The candidate spots can later be used as input of function
% "findSpotCentre1frame".
%
%
% INPUTS:
% frame: Input image, a single frame from a recorded sequence.
%
% method: two different methods, either 1 or 2.
%
% gaussian: a switch which runs a guassian filter if set to 1
%
% disk radius in pixels for top hat transformation, quite critical to the method for finding candidate bright spots.
%
%
% output=1 to see figures and advance manually
%
%
%
% OUTPUTS:
%result is matrix same size as the image with 1s at candidate spots
%

if output==1
    pause on
end

frame0 = mat2gray(frame);
if output==1
    imshow(frame0,[])
    title('original image frame')
    pause
end

if gaussian==1
    frame0=imfilter(frame0,fspecial('gaussian',3));
    if output==1
        imshow(frame0,[])
        title('gaussian filtered image')
        pause
    end
else
end

se = strel('disk',disk_radius); % structural element, disk of radius 5 pixels (default). The radius is quite critical to the method.
Signal = imtophat(frame0,se); % enhanced signal-image via "top hat" transformation (evens out the background).
if output==1
    imshow(Signal,[])
    title('enhanced signal-image via "top hat" transformation (evens out the background)')
    pause
end
% "Top-hat": substract an opened image from the original. Opening with a
% large enough structural element (disk of radius at least 5 pixels) ensures that
% the opening of the image produces a reasonable estimate of the background
% Bgnd0 = imopen(frame0,se); % Background estimate. Morphological opening (=erosion followed by dilation) with large enough structural element ("se" does not fit entirely within the foreground).
% if output==1
%     imshow(Bgnd0,[])
%     title('Background estimate. Morphological opening (=erosion followed by dilation) with large enough structural element ("se" does not fit entirely within the foreground).')
%     pause
% end
% % Middle row in image that separates top half (red channel) from bottom
% % half (green channel). If the image is 512x512, half = 256.
% half = round(size(frame0,1)/2); % round to closest integer.
% whole = round(size(frame0,1));
% % Separating foreground (signal) from background:
% 
% SignalMask = im2bw(Bgnd0,graythresh(Bgnd0)); % graythresh gives the threshold and im2bw converts image into a thresholded black and white (bw) image.
% 
% if output==1
%     imshow(SignalMask,[])
%     title('thresholded image')
%     pause
% end
% Signal3 = Signal.*SignalMask; % enhanced signal-only image.
% if output==1
%     imshow(Signal3,[])
%     title('enhanced signal-only image')
%     pause
% end
% Thresholding the signal-only image to try and find spot candidates:
%sigBW = im2bw(Signal3,graythresh(Signal3)); 
%sigBW = im2bw(Signal,graythresh(Signal)); 

[histy,histx]=imhist(Signal(Signal>0));

%hist(Signal(Signal>0))

[widthx, maxminvalue] = fwhm(histx,histy);
if widthx>0
threshold_cell_auto=maxminvalue+0.8*widthx; %this seems to kinda work for Adam Mig1-GFP data for cell outlines
threshold_cell_auto=maxminvalue;
else
   % threshold_cell_auto=maxminvalue;
    threshold_cell_auto=histx(histy==max(histy));
end
h = fspecial('gaussian');
GaussSignal=imfilter(Signal,h);
%GaussSignal2 = imtophat(GaussSignal,se);
GaussSignal2=GaussSignal;
sigBW = im2bw(GaussSignal2,threshold_cell_auto);
%sigBW = im2bw(GaussSignal2,graythresh(GaussSignal2));
if output==1
   % bar(histx, histy)
   % title(strcat('threshold=',num2str(threshold_cell_auto)))
   imshow(GaussSignal,[])
   title('Gaussian filtered tophat')
    pause
%       imshow(GaussSignal2,[])
%    title('Top-hat Gaussian filtered tophat')
%     pause
    imshow(sigBW,[])
  %  title('thresholded signal-only image')
  title('thresholded top-hat transformed image')
    pause
end

% filling holes here sort of works but leads to too many candidates
%sigBW=bwmorph(sigBW,'fill',1); 

B1 = imopen(sigBW,strel('disk',1)); % morph. opening with disk of radius 1 pixel.
B2 = bwmorph(B1,'fill',1); % fill holes of size one pixel, apply operation one time (1).
%B3 = bwmorph(B2,'open',1); % open image, apply operation one time.
B3=B2;
if output==1
%     imshow(sigBW,[])
%     title('fill holes of size one pixel - new step')
%     pause
    imshow(B1,[])
    title('morph. opening with disk of radius 1 pixel.')
    pause
    imshow(B2,[])
    title('fill holes of size one pixel, apply operation one time (1).')
    pause
%     imshow(B3,[])
%     title('open image, apply operation one time.')
%     pause
end

% Method 1:
if method==1
    B4 = imerode(B3, [1 1; 1 1]); % erode with a 2x2 square shape.
    B5 = imerode(B4, [1 1; 1 1]); % erode with a 2x2 square shape.
    B6 = bwmorph(B5,'shrink',1); % shrinks objects to points, apply operation one time.
    % B6 is a binary image (matrix) with ones at positions of candidate spots and zeros elsewhere.
    result = B6;
    if output==1
        imshow(B4,[])
        title('erode with a 2x2 square shape.')
        pause
        imshow(B5,[])
        title('erode with a 2x2 square shape. again')
        pause
        imshow(B6,[])
        title('B6 is a binary image (matrix) with ones at positions of candidate spots and zeros elsewhere.')
        pause
    end
    
end


% Method 2:
if method==2
    if max(B3(1:end))>0
    % B7 = bwdist(~B3); % distance transform (distance between each pixel and its nearest nonzero pixel) of the complement of the image
    B8 = bwulterode(B3); % "ultimate erosion" of a binary image = regional maxima of the Euclidean distance transform of the complement of the image.
    else
        B8=B3;
    end
    result = B8; % B8 is a binary image (matrix) with ones at positions of candidate spots and zeros elsewhere.
    if output==1
        imshow(B8,[])
        title('"ultimate erosion" of a binary image = regional maxima of the Euclidean distance transform of the complement of the image.')
        pause
    end
end
%result=B3;
if exist('CellMask')
    result=result.*SignalMask;
else
end


end