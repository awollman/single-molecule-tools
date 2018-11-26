function spotsBaseline=overTracker(SpotsCh1, imageData, noBaseFrames,p)

subarray_halfwidth = p.subarray_halfwidth; % (Default: 8 pixels). Halfwidth of image square subarray
% ROI, typically a square of size 17x17 pixels.
inner_circle_radius = p.inner_circle_radius; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
gauss_mask_sigma = p.gauss_mask_sigma; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.

SpotLastFrame=zeros(max(SpotsCh1(:,10)),1);
LFindex=zeros(max(SpotsCh1(:,10)),1);
y_estimate=zeros(max(SpotsCh1(:,10)),1);
x_estimate=zeros(max(SpotsCh1(:,10)),1);
spot_num=size(SpotsCh1,1)+1;
spotsBaseline=SpotsCh1;
%Loop over trajectory numbers
for k=1:max(SpotsCh1(:,10))
    [SpotLastFrame(k), LFindex(k)]=max(SpotsCh1(SpotsCh1(:,10)==k,9));
    SpotsCh1X=SpotsCh1(SpotsCh1(:,10)==k,1);
    SpotsCh1Y=SpotsCh1(SpotsCh1(:,10)==k,2);
    y_estimate(k,1)=SpotsCh1Y(end);
    x_estimate(k,1)=SpotsCh1X(end);
    %Loop over frames
    if SpotLastFrame(k)+noBaseFrames < max(SpotsCh1(:,9))
    for p=SpotLastFrame(k)+1:SpotLastFrame(k)+1+noBaseFrames
        frame=imageData(:,:,p);
        [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
            findSpotCentre2noloopLocal(frame,x_estimate(k),y_estimate(k),subarray_halfwidth,inner_circle_radius,gauss_mask_sigma,0.05, 0,0);
        snr1=Isp/(bg_noise_std*mask_pixels);
        spotsBaseline(spot_num,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, inner_circle_radius, inner_circle_radius, -1, p,k, snr1, SpotsCh1(1,12)];
        spot_num=spot_num+1;
    end
    end
    
   


    
end
end

 function [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std,  mask_pixels, noConvergenceFlag]= ...
    findSpotCentre2noloopLocal(imageData,x_estimate,y_estimate,subarray_halfwidth,inner_radius,sigma_gauss, error_set, clip_override,show_output)

% make imageData of class double:
imageData = double(imageData);

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
if (round(y_estimate)+d)>size(imageData,1)
    y_estimate = size(imageData,1)-d;
end
if (round(x_estimate)-d)<1
    x_estimate = d+1;
end
if (round(x_estimate)+d)>size(imageData,2)
    x_estimate = size(imageData,2)-d;
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
I = imageData(round(y_estimate)-d:round(y_estimate)+d,round(x_estimate)-d:round(x_estimate)+d);
Idata=I;


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
    

    % The background mask is just the negative of the inner circle mask, i.e.,
    % not-inner_mask_0. It has zeros in centre and ones in surrounding pixels.
    bgnd_mask = double(~inner_mask); % make of class double (not binary) for what follows.
    
    pos_bgnd = find(bgnd_mask==1); % positions of bgnd only intensities in bgnd mask.
    Ibgnd = I.*bgnd_mask; % bgnd region.
    
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
      
    noConvergenceFlag=1;
    

end

