%Threshold's images using one of three methods:
%method 0, uses Otsu's method
%method 1, uses fixed threshold specified in input
%method 2, uses the FWHM of image histogram
%method 3, threshold is a real pixel value
%method 4, otsus with an adaptive threshold

%% Ouputs
%Mask is a stack of binary masks
%% Inputs
%image=grayscale image
% sizeThresh= size threshold below which objects are removed
%threshold=actual threshold if method 3
%%
function [Mask]=thresholdSegment(image,sizeThresh,threshold_method,threshold)
  %  
  
  switch threshold_method
      case 1
          image=mat2gray(image);
      threshold=threshold;
      bw = im2bw(image, threshold); 
      case 2
          image=mat2gray(image);
%           binNum=max(image(1:end))-min(image(1:end));
%    [histy,histx]=hist(double(image(1:end)),double(binNum));
     [histy,histx]=imhist(image);
       [widthx, maxminvalue] = fwhm(histx,histy);
        threshold=(maxminvalue+widthx)/max(histx);
        bw = im2bw(image, threshold); 
      case 0
          image=mat2gray(image);
      threshold=graythresh(image);
      if threshold==0
          disp('oops, Otsus method failed big time')
      end
      bw = im2bw(image, threshold); 
      case 3
          bw = image>threshold;
      case 4
          bw=imbinarize(image,'adaptive');
  end
    
    %clean that up and then overlay the perimeter on the original image.
    bw2 = imfill(bw,'holes'); %fills in holes in image by flooding
    se = strel('disk',2);%disk se
    bw3 = imopen(bw2, se );%morphological opening to remove objects smaller than se 
    bw4 = bwareaopen(bw3, sizeThresh); %Remove small objects from binary image here with <50pixels total (inner ROI has 5 pixel raidus, so take 4 whole pixels radius as limiting size)
   
cc = bwconncomp(bw4,4);
L_cells = labelmatrix(cc);
% Mask=L_cells;
    Mask=zeros(size(image,1),size(image,2),max(L_cells(:)));
   % clear Mask
    for lbl = 1:max(L_cells(:)); %cycle through all objects in labe matrix... asumume only one nucleus per cell... could fall down here if problem with extended maxima of course!   

            Mask(:,:,lbl) = L_cells == lbl; %# find pixels belonging to current label. First object is the borders I think, so take the objects after these as indivudal cell masks  
    end
    Mask=Mask;
    

end