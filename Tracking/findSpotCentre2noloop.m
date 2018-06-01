function [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std,  mask_pixels, noConvergenceFlag]= ...
    findSpotCentre2noloop(image_data,x_estimate,y_estimate,subarray_halfwidth,inner_radius,sigma_gauss, error_set, clip_override,show_output)
% 
% Will return initial values if NO CONVERGE
%
% Created by Isabel Llorente-Garcia, (2011) with a little help from Adam
% Wollman in 2014
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Function to iteratively find the centre of a fluorescence spot on an image.
% 
% 
% Inputs:  
% image_data is a matrix containing the image data for a single frame. 
% image_data is a matrix with original values, so not the grayscale version between 0 and 1.
% image_data is later made of class double.
% It can be the output of extract1frame(frame_number) or extract_image_sequence_data(image_label), for example.
%
% x_estimate, y_estimate are the initial estimated centre positions for the
% iterative method.
%
% subarray_halfsize: the algorithm selects this number of pixels above and below and to
% left and right of estimated centroid pixel to form the selected image
% subarray (I). (The default value is 8 pixels).
%
% inner_radius: radius of inner circular mask in pixels (default is 5 pixels). This mask moves
% around inside the subarray to find the bright spot.
%
% sigma_gauss: sigma of gaussian mask in pixels (default is 2 pixels).
%
% guess_sigma_fit: starting guess for Gaussian fit of bright
% spot intensity (default = 3).
%
%clip_override=1 if you want to return initial data if clipping flag=1
%
% Outputs: values resulting from the iterative method. The output is a
% structure "spot" with the following fields (see end of file).
%
% Example of how to call this function: s1 = findSpotCentre1frame(A,271,364,8,5,2), 
% where A is the image array (single frame data). Then to obtain results
% from the structure s1, do, for example: s1.CentreX, s1.CentreY,
% s1.SigmaFit, s1.I0Fit, etc.

% make image_data of class double:
image_data = double(image_data);

d = subarray_halfwidth; % subarray halfsize. d needs to be at least equal to the radius of the inner mask, inner_radius. Default: d=8, inner_radius=5.
if show_output==1
    pause on
end
%-----------------------
% error control:
if subarray_halfwidth < inner_radius
    error('findSpotCentre1frame:one','check input parameters: "subarray_halfsize" cannot be smaller than "inner_radius"')
end

% If the spot candidate is at the edge of the image move it away from the edge until we can take a
% subarray around it:
if (round(y_estimate)-d)<1
    y_estimate = d+1;
end
if (round(y_estimate)+d)>size(image_data,1)
    y_estimate = size(image_data,1)-d;
end
if (round(x_estimate)-d)<1
    x_estimate = d+1;
end
if (round(x_estimate)+d)>size(image_data,2)
    x_estimate = size(image_data,2)-d;
end    
%------------------------

% PARAMETERS FOR FINAL GAUSSIAN FITTING OF SPOT:
% Set fit constrains: upper and lower bounds for fit parameters I0 and sigma_fit:
min_I0 = 0; 
max_I0 = Inf; % unconstrained.
min_sigma_fit = -inner_radius; % rather generous limits here.
max_sigma_fit = inner_radius; % rather generous limits here.

% -----------------------

% create image subarray (I) of size (2*d+1)x(2*d+1) (default size 17x17) centred on the
% centroid estimate pixel (x_estimate,y_estimate). I is fixed during the iterative process of finding the centre of the bright spot.
I = image_data(round(y_estimate)-d:round(y_estimate)+d,round(x_estimate)-d:round(x_estimate)+d);
Idata=I;
if show_output==1
    imshow(I, [],'InitialMagnification',1000)
    title('image subarray (I) of size (2*d+1)x(2*d+1) (default size 17x17) centred on thecentroid estimate pixel (x_estimate,y_estimate).')
    pause
end

% Create matrices containing the x and y positions corresponding to the
% subarray I. Xs is a matrix of the same size as I, containing x values for
% all pixels and similarly for Ys.
[Xs,Ys] = meshgrid(round(x_estimate)-d:round(x_estimate)+d,round(y_estimate)-d:round(y_estimate)+d);


% Initialisation values:
k = 1; % index for iteration loop to find centre position of bright spot.
error1 = 1; % initialisation value to get inside loop.
clipping_flag = 0;  % initialise to zero to get inside loop.
x_centre = x_estimate;
y_centre = y_estimate;
list_Xcentre = [x_estimate];
list_Ycentre = [y_estimate];
list_deltaX = [];
list_deltaY = [];
list_error1 = [];
list_iteration = [0];
list_clip_flag = [0];
list_Ibg = []; % Average background intensity in image subarray.
list_Isp = []; % Total intensity within the spot (inner circle) mask in subarray, background corrected.
list_Itot = []; % Total subarray image intensity, background corrected. 
% The previous lists are there to keep track of the convergence of the method over
% the several iterations.


    %
    % % Size of mask matrix is 2*radius+1. All elements in fspecial('disk',5) are between 0 and 1 and the sum of all elements is 1.
    % % the ~=0 gives matrix with 1 at positions of elements different from 0, so it gives a circular mask of ones.
    % % inner_mask_0 = fspecial('disk',5)~=0;
    % % inner_mask_0 = padarray(inner_mask_0,[3,3]); % pad the previous matrix
    % % with zeros around it to get to a matrix of the same size as I.

    % Create inner circle mask: circle of radius 5 pixels, with ones inside and zeros outside.
    % This mask can move around inside the fixed image subarray I.
    % In the same loop (to save time): Create Gaussian mask of sigma=3pixels, centred on (x_centre,y_centre),
    % with matrix size (2*d+1)x(2*d+1) (default size 17x17).                    ((mask_gauss = fspecial('gaussian',17,3);))
    inner_mask = zeros(size(I,1)); % pre-define inner_mask size.
    % radius of inner circle mask is inner-radius. Default is 5 pixels.
    mask_gauss = zeros(size(I,1)); % pre-define mask_gauss size.
    % default fixed gaussian width is sigma_gauss = 3 pixels. 
    for ii = 1:size(I,1) % ii = index for rows
        for jj = 1:size(I,2) % jj = index for columns
            % inner circle mask: if distance to point (x_centre,y_centre) is <=inner_radius, then 1, otherwise, 0:
            sum_sq = (Xs(ii,jj)-x_centre)^2 + (Ys(ii,jj)-y_centre)^2;
            if round(sqrt(sum_sq))<=inner_radius 
                inner_mask(ii,jj)=1; 
            else
                inner_mask(ii,jj)=0; 
            end
            % gaussian mask:
            mask_gauss(ii,jj) = exp(-sum_sq/(2*sigma_gauss^2)); 
        end
    end
    mask_pixels=sum(sum(inner_mask));
    mask_gauss = mask_gauss/sum(sum(mask_gauss)); % normalise the Gaussian mask at the end.
    if show_output==1
        imshow(mask_gauss, [],'InitialMagnification',1000)
        title('Gaussian Mask')
    end
    

    % The background mask is just the negative of the inner circle mask, i.e.,
    % not-inner_mask_0. It has zeros in centre and ones in surrounding pixels.
    bgnd_mask = double(~inner_mask); % make of class double (not binary) for what follows.
    
    pos_bgnd = find(bgnd_mask==1); % positions of bgnd only intensities in bgnd mask.
    Ibgnd = I.*bgnd_mask; % bgnd region.
    if show_output==1
        imshow(Ibgnd, [],'InitialMagnification',1000)
        title('bgnd region.')
        pause
    end
    
    % ----------------------------------------
    % Total background intensity in background region before background correction:
    % Ibg_tot = sum(sum(I.*bgnd_mask)); % this gives same result.
    Ibg_tot = sum(Ibgnd(pos_bgnd));
    % Median background intensity in bgnd region, per pixel, Ibg_avg.
    % The number of pixels in the background mask equal to 1 is
    % length(find(bgnd_mask~=0)):
    % Ibg_avg = Ibg_tot/length(find(bgnd_mask~=0)); % this is a scalar.
%     Ibg_avg = median(Ibgnd(pos_bgnd)); % use median instead of mean for the bgnd to exclude hot pixels or nearby bright spots.
    Ibg_avg = mean(Ibgnd(pos_bgnd)); % use mean for the bgnd to exclude hot pixels or nearby bright spots.
    % Total intensity in inner mask region before background correction:
    Iinner_tot = sum(sum(I.*inner_mask));
    %
    % Calculate background-corrected subarray image:
    I2 = I-Ibg_avg;
    if show_output==1
        imshow(I2, [],'InitialMagnification',1000)
        title('background-corrected subarray image:')
        pause
    end
    
    % Calculate standard deviation of remaining noise in background
    % corrected image I2:
    bg_noise_offset_afterBGsubtract = mean(I2(pos_bgnd)); % offset noise level per pixel after background subtraction, should be close to zero.
    bg_noise_std = std(I2(pos_bgnd)); % standard deviation of matrix elements in bgnd region.
    % Total spot intensity (within the inner circle mask), background corrected:
    Isp = sum(sum(I2.*inner_mask));     
    % -----------------------------------------
    
    % Calculate product of background-corrected image and Gaussian mask:
    I3 = I2.*mask_gauss;
        if show_output==1
        imshow(I3, [],'InitialMagnification',1000)
        title('product of background-corrected image and Gaussian mask')
        pause
        end
    noConvergenceFlag=1;
    
%     % Calculate revised estimates of x and y centre positions by multiplying I2 by
%     % a Gaussian mask (result is I3) and weighing positions by the intensity value in that result:
%     x_centre_new = sum(sum(I3.*Xs))/sum(sum(I3));
%     y_centre_new = sum(sum(I3.*Ys))/sum(sum(I3));
%     
%     delta_x = x_centre_new - x_centre; % rounded to integer pixels
%     delta_y = y_centre_new - y_centre;
%     
%     % error1 is the error deviation used to decide if the method has converged or not:
%     error1 = sqrt(delta_x^2 + delta_y^2); % error 1 is the distance in pix between old and new centroid estimates.
%     % error1 = abs(delta_x) + abs(delta_y); % this is an alternative error that could be considered.
%     
%     % disp(['delta x = ',num2str(delta_x)])
%     % disp(['delta y = ',num2str(delta_y)])
%     % disp(['error deviation = ',num2str(error1)])
%     
%     % -----------------------------------------------------
%     % Clipping flag: 
%     % If the inner circle mask moves too far away from the fixed centre of
%     % the image subarray, i.e., if the inner-circle mask clips the edge of
%     % the square subarray, the "Clipping flag" goes on and takes a value of
%     % 1. Otherwise, it is zero. A clipping flag value of 1 indicates the
%     % found spot is not very reliable.
%     d_found_centre_x = abs(x_centre_new - x_estimate); 
%     d_found_centre_y = abs(y_centre_new - y_estimate); 
%     
%     if d_found_centre_x >(subarray_halfwidth-inner_radius+1) || d_found_centre_y >(subarray_halfwidth-inner_radius+1) % inner mask clips edge of subarray.
%         clipping_flag = 1;
%     else
%         clipping_flag = 0;
%     end
%     % -----------------------------------------------------
%     
%     % update lists:
%     list_Xcentre(k+1) = x_centre_new;
%     list_Ycentre(k+1) = y_centre_new;
%     list_deltaX(k) = delta_x;
%     list_deltaY(k) = delta_y;
%     list_error1(k) = error1;
%     list_clip_flag(k) = clipping_flag;
%     list_iteration(k+1) = k;
%  
%     list_Ibg(k) = Ibg_avg; % Average background intensity in image subarray.
%     list_Isp(k) = Isp; % Total intensity within the spot (inner circle) mask in subarray, background corrected.
%     if k==1
%         Isp1=Isp;
%          Ibg_avg1= Ibg_avg;
%     end
%     % re-asign values for next iteration:
%     x_centre = x_centre_new;
%     y_centre = y_centre_new;
%     k = k+1;
%     
%     % prevent infinite loop:
%     if k>300
%         % disp('The method did not reach convergence after 300 iterations.')
%         noConvergenceFlag = 1;
%         % Note that in this case the following break gets us out of the
%         % loop.
% 
%     else
%         noConvergenceFlag = 0;
%     end   
    
% 
% if noConvergenceFlag == 1
%     x_centre=x_estimate;
%     y_centre=y_estimate;
% %     clipping_flag=0;
%     Ibg_avg=Ibg_avg1;
%     Isp=Isp1;
% else
%      x_centre = x_centre_new;
%     y_centre = y_centre_new;
% end
% if clipping_flag==1
%     if clip_override==1
%         x_centre=x_estimate;
%         y_centre=y_estimate;
%         %     clipping_flag=0;
%         Ibg_avg=Ibg_avg1;
%         Isp=Isp1;
%     end
% end
end

