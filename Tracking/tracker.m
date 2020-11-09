function [SpotsCh1, SpotsCh2, frame_average,p, meta_data, image_data,spotImages] = tracker(fileName,p)

% Code for tracking bright foci in image stacks
%INPUTS
%filename: name of image stack or folder containing series of tifs
%p: parameter structure (see below for details) if not set default values
%used
%OUTPUTs
% SpotsCh1/2 is an array, each row contains the information for a spot found in an image frame in the series. The columns contain the following information:
% 
% 1. X coordinate (pixels)
% 2. Y coordinate (pixels)
% 3. Clipping_flag (a switch, please ignore)
% 4. Mean local background pixel intensity
% 5. Total spot intensity, background corrected
% 6. X spot sigma width
% 7. Y spot sigma width
% 8. Spot central intensity (peak intensity in a fitted Gaussian)
% 9. Frame number the spot was found in
% 10. Trajectory number, spots in the same trajectory have the same trajectory number
% 11. Signal to noise ratio
% 12. Frame in which laser exposure began
% 
%frame_average: frame average image
%p: output parameters
% meta_data: any meta_data stored in the image
%image_data: the raw imaging data tracked
% spotImages: every spot tracked - NOT RECOMMENDED TO USE

if exist('p','var')==1
    %Read in parameters
else
    % Number of frames to average over when calculating a frame average in
    % used  as image to choose spots in cursor
    % mode and output
p.noFrames=5; % number of frames to average over, starting from start_frame. (default: 5, end_frame-start_frame)
    % PARAMETERS for finding spot centres
    % The total integrated spot intensity is a bgnd corrected one, inside a
    % circular mask of radius inner_circle_radius.
    p.subarray_halfwidth = 8; % (Default: 8 pixels). Halfwidth of image square subarray
    %p.subarray_halfwidth = 7; % (Default: 8 pixels). Halfwidth of image square subarray
    % ROI, typically a square of size 17x17 pixels.
    p.inner_circle_radius = 5; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
    %p.inner_circle_radius = 3; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
    p.gauss_mask_sigma = 2; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.
    p.guess_sigma_Fit = 3; % starting guess for Gaussian fit of brightspot intensity (default = 3).
    
    % PARAMETERS for deciding if we accept a spot centre
    p.sigmaFit_min = 0;  % minimum acceptable sigma of gaussian fit to spot, in pixels (2) (-3).
    p.sigmaFit_max = p.inner_circle_radius; % maximum acceptable sigma of gaussian fit to spot, in pixels (4) (3).
    p.SNR_min = 0.4; % minimum acceptable signal-to-noise ratio (at least 0.4)
    
    % PARAMETERS for eliminating coincident spots:
    p.d_min= 1; % distance (in pixels) for eliminating coincidences
    
    % PARAMETERS for building trajectories:
    % For linking spots in current and previous frames:
    p.d_01_max = 5; % max distance in pixels between spot centres in current and previous frames, for linking them into a trajectory (5).
    p.Iratio_01_min = 0.5; % min ratio of total spot intensities (after bgnd subtraction) (0.5).
    p.Iratio_01_max = 3; % max ratio of total spot intensities (after bgnd subtraction) (frame k-1/frame k) (large enough value (3) to account for blinking).
    p.SigmaRatio_01_min = 0.5; % min ratio of spot widths (sigma of Gaussian fit) (0.5).
    p.SigmaRatio_01_max = 3; % max ratio of spot width (sigma of Gaussian fit) (2).
    p.error_set=0.05; %error in iterative Gaussian masking
    p.exclude_region = 2; % Parameter to exclude a region from the middle if Csplit>0;
    p.start_channel=1; %Channels to use
    p.end_channel=1;
    p.disk_radius = 5; % for finding spots in image
    p.ALEX=0; % Switch, if ALEX experiment=1
    
    p.satPixelVal=10^10; % Checks for pixel values above this value and won't track, currently set to crazy value, =15898 on Andor 128
    
    
    % There are 2 methods for finding candidate spots which work better for
    % different datasets or if set =3 runs both and uses all spots
    % (recommended)
    p.CandidateFindMethod=3;
    
    % Choose to iterate a Gausian to determine PSF width (1D only) ==1
    % or explicitly fit a 2D Gaussian using MAtlab fitting routines==2
    % or fit a fully rotating Gaussian==3
    p.GaussSwitch=1;
    % switch to use parallel processing or not
    p.useParallel=0;
    %Distance to remove candidates which are too close together, set to 0 to
    %ignore
    p.Candidate_d_min=0;
    
    %set to one and saves an array with every spot image saved
    p.spotImageSave=0;
    %gaussian=1 if running a gaussian filter over image data before finding
    %spots
    p.gaussian=0;
    
    % Use this to specify interesting spots with the cursor, rather than
    % autodetecting them
    p.use_cursor=0;
    
    % Use this to open files with BioFormats to open non tiff files
    p.useBioFormats=0;
    p.all=1; %Use this keyword to load entire image file
    p.startFrame=1; % Or specify start and end frames
    p.endFrame=999;
    p.firstLeft=2; %If determine first frames=0 these are used
    p.firstRight=1;
    
    
    %Specify how many frames to track after the laser has switched on
    p.FramesToTrack=0;
    
    % Set this to 1 if there are blank frames before laser turns on or shutter
    % opens
    p.DetermineFirstFrames=0;
    
    % Switch, =1 to determine laser on time with differential rather than max
    % intensity
    p.use_diff=0;
    % CSplit defines how the channels are split, =0 for whole frame (no split),
    % 1 for left/right and 2 for up/down
    p.CSplit=0;
    %If this =1, then graphs will appear
    p.show_output=0;
    % set this and graphs will appear for every spot!
    p.show_all_output=0;
    p.show_text_output=1;
    
end
%Initialise Spot variables
SpotsCh1=[];
SpotsCh2=[];

%% CREATE FITTYPE AND OPTIONS
switch p.GaussSwitch
    case 1
        % don't need any fitting stuff
        myfit=[];
        options=[];
    case 2
% Create fit type for constrained 2D Gaussian fit to spots
myfit = fittype('Ibg_avg+(Isp./(2.*pi.*sdx.*sdy))*exp(-(((x-x0).^2)./(2.*sdx^2)+((y-y0).^2)./(2.*sdy^2)))',...
    'problem', {'Ibg_avg','Isp','x0','y0'}, 'independent', {'x', 'y'}, 'dependent', 'z');
% Fit options:
options = fitoptions(myfit);
options.StartPoint = [p.guess_sigma_Fit, p.guess_sigma_Fit];
options.Lower = [p.sigmaFit_min, p.sigmaFit_min];
options.upper = [p.sigmaFit_max, p.sigmaFit_max];
    case 3
myfit = fittype('Ibg_avg+ (Isp./(2.*pi.*sdx.*sdy))*exp( - ((cos(theta)^2/(2*sdx^2) + sin(theta)^2/(2*sdy^2))*(x-x0).^2 - 2*(-sin(2*theta)/(4*sdx^2) + sin(2*theta)/(4*sdy^2))*(x-x0).*(y-y0) + (sin(theta)^2/(2*sdx^2) + cos(theta)^2/(2*sdy^2))*(y-y0).^2))',...
    'problem', {'Ibg_avg','Isp','x0','y0'}, 'independent', {'x', 'y'}, 'dependent', 'z');
options = fitoptions(myfit);
options.StartPoint = [p.guess_sigma_Fit, p.guess_sigma_Fit,0];
options.Lower = [p.sigmaFit_min, p.sigmaFit_min,0];
options.upper = [p.sigmaFit_max, p.sigmaFit_max,2*pi];
end


%% OPEN DATA

%Open tif data from image_label which is whatever string is before the file
%extention
if p.useBioFormats==1
    if isa(fileName,'char')==1
        if p.all==1
            [image_data, meta_data]=imEx1(fileName);
        else
            [image_data, meta_data]=imEx1(fileName,p.startFrame,p.endFrame);
        end
        numFrames=size(image_data,3);
        
    else
        image_data=fileName;
        numFrames=size(image_data,3);
    end
else
    if isa(fileName,'char')==1
        [numFrames, ~, ~, image_data, ~] = ExtractImageSequence3(fileName, p.all, p.startFrame, p.endFrame);
        meta_data='sorry metadata not available in this mode';
    else
        image_data=fileName;
        numFrames=size(image_data,3);
        meta_data='sorry metadata not available in this mode';
    end
end
%[numFrames, frame_Ysize, frame_Xsize, image_data, image_path] = extract_image_sequence_dataAWarray(image_label, all);
disp('data loaded')

%% DETERMINE LASER ON FRAME

%Determine when laser turned on, detects both channels seperately if ALEX
if p.DetermineFirstFrames==1
    [firstLeft, firstRight, ~, ~] = LaserOn3(image_data, p.use_diff,p.ALEX);
    disp('start determined')
else
    firstLeft=p.firstLeft;
    firstRight=p.firstRight;

end

%Number of frames to track
%FramesToTrack=numFrames-firstLeft-1;
%FramesToTrack=50;
%firstLeft=firstRight;
%firstLeft=1;

if p.FramesToTrack==0
    p.FramesToTrack=numFrames-firstLeft;
end

%% CALCULATE FRAME AVERAGE

%Calculate Frame Average of the data
try
    if firstLeft<(size(image_data,3)-5)
        frame_average = FrameAverage3(image_data, p.noFrames, firstLeft,firstRight,p.ALEX);
        disp('average calculated')
    else
        disp('WARNING FRAME AVERAGE NOT CALCULATED AS STARTFRAME TOO CLOSE TO ENDFRAME')
        frame_average=image_data(:,:,firstLeft);
    end
catch err
    disp('calculating frame average failed because:')
    disp(err.message)
end

%%
for Ch=p.start_channel:p.end_channel
    % Initialise spot array
    spots=[];
    spotImages=[];
    
    
    if p.ALEX==0
        if Ch==1
            startFrame=firstLeft;
            endFrame=firstLeft+p.FramesToTrack;
        else
            startFrame=firstRight;
            endFrame=firstRight+p.FramesToTrack;
        end
        FrameInt=1;
    else
        FrameInt=2;
        if Ch==1
            startFrame=firstLeft;
            endFrame=firstLeft+p.FramesToTrack;
        else
            startFrame=firstRight;
            endFrame=firstRight+p.FramesToTrack;
        end
    end
    if endFrame>size(image_data,3)
        endFrame=size(image_data,3);
    end
    
    %Loop over frames
    for i=startFrame:FrameInt:endFrame
        if p.show_text_output==0
            %    h=waitbar(i-startFrame/p.FramesToTrack);
        end
        disp(strcat('Tracking frame: ',num2str(i)))
        %% FIND CANDIDATE SPOTS
        % Divide image into two channels, left and right
        if Ch==1
            if p.show_text_output==1
                disp('Ch1')
            end
            switch p.CSplit
                case 0
                    frame=image_data(:,:,i);
                case 1
                    frame=image_data(:,1:round(size(image_data,2)/2-p.exclude_region),i);
                case 2
                    frame=image_data(1:round(size(image_data,1)/2-p.exclude_region),:,i);
            end
            SpotsCh1=[];
        elseif Ch==2
            if p.show_text_output==1
                disp('Ch2')
            end
            switch p.CSplit
                case 0
                    frame=image_data(:,:,i);
                case 1
                    frame=image_data(:,round(size(image_data,2)/2+p.exclude_region):end,i);
                case 2
                    frame=image_data(round(size(image_data,1)/2-p.exclude_region):end,:,i);
            end
            SpotsCh2=[];
        end
        %% Check for saturation
        if max(max(frame))>p.satPixelVal
            continue
        end
        
        % For cursor mode
        if p.use_cursor==1
            if i==startFrame
                pause on
                % Display frame average to choose spots
                if Ch==1
                    frame_averageCH=frame_average(:,1:round(size(image_data,2)/(p.end_channel-p.start_channel+1)-p.exclude_region));
                    %  frame_averageCH=frame_average;
                elseif Ch==2
                    frame_averageCH=frame_average(:,round(size(image_data,2)/2+p.exclude_region):end);
                end
                h1=figure;
                imshow(frame_averageCH,[],'InitialMagnification','fit')%HM modify magnification
                title('click a spot and hold alt key to select multiple spots, push any key when finished')
                datacursormode on
                pause
                dcm_obj = datacursormode(h1);
                info_struct = getCursorInfo(dcm_obj);
                %Loop over spots chosen and pull out co-ordinates
                for q=1:size(info_struct,2)
                    Spot_coords=info_struct(q).Position;
                    y_estimate(q,1)=Spot_coords(2);
                    x_estimate(q,1)=Spot_coords(1);
                end
                close all
            end
        else
            % Create matrix of 1s where spots might be
            % Now 3 methods for doing this
            switch p.CandidateFindMethod
                case 1
                    [result] = findSpots2(frame,2,p.disk_radius,p.gaussian,0);
                case 2
                    [result] = findSpots3(frame,2,p.disk_radius,p.gaussian,0);
                case 3
                    [result1] = findSpots2(frame,2,p.disk_radius,p.gaussian,0);
                    [result2] = findSpots3(frame,2,p.disk_radius,p.gaussian,0);
                    result=result1+result2;
                    result(result>1)=1;
            end
            % Convert those to spot co-ordinates
            [y_estimateTemp, x_estimateTemp]=ind2sub(size(result), find(result));
            %  [y_estimate, x_estimate]=ind2sub(size(result), find(result));
            % Remove candidates which are too close together to yield
            % seperate spots-
            if p.Candidate_d_min>0
                [y_estimate,x_estimate]=MergeCoincidentCandidates2(y_estimateTemp, x_estimateTemp, p.Candidate_d_min);
            else
                y_estimate=y_estimateTemp;
                x_estimate=x_estimateTemp;
            end
            
        end
        if p.show_text_output==1
            disp('candidates found')
        end
        %Plot the candidate spots
        if p.show_output==1
            imshow(frame, [],'InitialMagnification','fit')%HM modify magnification
            hold on
            plot(x_estimate, y_estimate, 'o')
            hold off
            title('candidate spots')
            pause
        end
        %% FIT TO SPOTS AND REJECT
        %Loop over all found spots
        spots_temp=zeros(size(x_estimate,1),12);
        spotImageTemp=zeros(p.subarray_halfwidth*2+1,p.subarray_halfwidth*2+1,size(x_estimate,1));
        if p.useParallel==1 % use a parfor loop
            parfor j=1:size(x_estimate,1)
                if p.use_cursor==1
                    trajNo=j;
                    clip_override=1;
                else
                    trajNo=0;
                    clip_override=0;
                end
                % Iterative gaussian masking to determine spot centre
                [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
                    gaussianMask(frame,x_estimate(j),y_estimate(j),p,clip_override);
                % Calculate spots signal to noise ratio
                snr1=Isp/(bg_noise_std*mask_pixels);
                
                % check if already spot with those co-ords here
                % isSpotAlready=sum(((spots_temp(:,1)-x_centre).^2+(spots_temp(:,2)-y_centre).^2).^0.5<p.d_min);
                isSpotAlready=0;
                % If in curosr mode OR a spot was found and it didn't clip
                if p.use_cursor==1 || noConvergenceFlag==0 && clipping_flag==0
                    % Only store spots with good enough snr
                    if p.use_cursor==1 || snr1>p.SNR_min && isSpotAlready<1
                        switch p.GaussSwitch
                            case 1
                                [sdx, sdy, Icent] = iterate1DgaussianFixedCenter3(frame,Ibg_avg, Isp, x_centre, y_centre,p);
                                spots_temp(j,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];
                            case 2
                                [sdx, sdy, Icent,Rsquare] = fit2DgaussianFixedCenter4(frame,Ibg_avg, Isp, x_centre, y_centre,myfit,options,p);
                                spots_temp(j,:)=[x_centre, y_centre, Rsquare, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];
                            case 3
                                [sdx, sdy, Icent,theta] = fit2DRotGaussianFixedCenter3(frame,Ibg_avg, Isp, x_centre, y_centre,myfit,options,p);
                                spots_temp(j,:)=[x_centre, y_centre, theta, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];
                            otherwise
                        end
                        % The spot array, 10th field is trajectory number,
                        % initialised to 0
                        spotImageTemp(:,:,j)=Idata;
                        
                    end
                end
                
            end
            spotImageTemp(:,:,spots_temp(:,1)==0)=[];
            spots_temp(spots_temp(:,1)==0,:)=[]; 
            % get rid of identical spots, within ~d_min of each other
            
            if isempty(spots_temp)==0
                [~,ia,~]=uniquetol(spots_temp(:,1)+spots_temp(:,2),p.d_min/max(spots_temp(:,1)+spots_temp(:,2)));
              %  [~,ia,~]=unique(spots_temp(:,1)+spots_temp(:,2));
       %       try
                spots_temp2=spots_temp(ia,:);
%               catch
%                   spots_temp(ia,:)
%                   spots_temp
%               end

                spotImageTemp2=spotImageTemp(:,:,ia);

            end
            
            % wipe the stored spot image array if spotImageSave=0 to save
            % memory
            if p.spotImageSave==0
                spotImageTemp2=[];
            end
            
            if exist('spots_temp2')
            spots=cat(1,spots,spots_temp2);
            spotImages=cat(3,spotImages,spotImageTemp2);
            end
        else %use a regular for loop
            for j=1:size(x_estimate,1)
                if p.use_cursor==1
                    trajNo=j;
                    clip_override=1;
                else
                    trajNo=0;
                    clip_override=0;
                end
                % Iterative gaussian masking to determine spot centre
                [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
                    gaussianMask(frame,x_estimate(j),y_estimate(j),p,clip_override);
                % Calculate spots signal to noise ratio
                snr1=Isp/(bg_noise_std*mask_pixels);
                
                % check if already spot with those co-ords here
                % isSpotAlready=sum(((spots_temp(:,1)-x_centre).^2+(spots_temp(:,2)-y_centre).^2).^0.5<p.d_min);
                isSpotAlready=0;
                % If in curosr mode OR a spot was found and it didn't clip
                if p.use_cursor==1 || noConvergenceFlag==0 && clipping_flag==0
                    % Only store spots with good enough snr
                    if p.use_cursor==1 || snr1>p.SNR_min && isSpotAlready<1
                        switch p.GaussSwitch
                            case 1
                                [sdx, sdy, Icent] = iterate1DgaussianFixedCenter3(frame,Ibg_avg, Isp, x_centre, y_centre,p);
                                                        spots_temp(j,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];

                            case 2
                                [sdx, sdy, Icent,Rsquare] = fit2DgaussianFixedCenter4(frame,Ibg_avg, Isp, x_centre, y_centre,myfit,options,p);
                                spots_temp(j,:)=[x_centre, y_centre, Rsquare, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];

                            case 3
                                [sdx, sdy, Icent,theta] = fit2DRotGaussianFixedCenter3(frame,Ibg_avg, Isp, x_centre, y_centre,myfit,options,p);
                                spots_temp(j,:)=[x_centre, y_centre, theta, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, startFrame];
                            otherwise
                        end
                        % The spot array, 10th field is trajectory number,
                        % initialised to 0
                        
                          spotImageTemp(:,:,j)=Idata;
                  %      spots_temp(j,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,trajNo, snr1, firstLeft];
                    end
                end
                
            end
            spotImageTemp(:,:,spots_temp(:,1)==0)=[];

            spots_temp(spots_temp(:,1)==0,:)=[]; %HM, gets rid of empty rows
            % get rid of identical spots, within ~d_min of each other
            
            if isempty(spots_temp)==0
                [~,ia,~]=uniquetol(spots_temp(:,1)+spots_temp(:,2),p.d_min/max(spots_temp(:,1)+spots_temp(:,2)));
                spots_temp2=spots_temp(ia,:);
                spotImageTemp2=spotImageTemp(:,:,ia);
            end
            
          if p.spotImageSave==0
                spotImageTemp2=[];
            end
            
            if exist('spots_temp2')
            spots=cat(1,spots,spots_temp2);
              spotImages=cat(3,spotImages,spotImageTemp2);
            end
        end
        
        
        if isempty(spots)==0
            if p.show_text_output==1
                disp('Centres found and fitted')
            end
            
            %% PLOT ALL FOUND SPOTS ON IMAGE
            %Plot all the found spots on image
            
            if p.show_output==1
                imshow(frame,[],'InitialMagnification','fit')%HM modify magnification
                hold on
                title('Found elipses on original image')
                %    plot(spots(spots(:,9)==i,1),spots(spots(:,9)==i,2), 'o')
                %Plot elipses on original image
                if p.GaussSwitch==3
                    h=ellipse(spots(spots(:,9)==i,6),spots(spots(:,9)==i,7),spots(spots(:,9)==i,3),spots(spots(:,9)==i,1),spots(spots(:,9)==i,2),'b');
                    for k=min(find(spots(:,9)==i)):max(find(spots(:,9)==i))
                        text(spots(k,1)+3,spots(k,2)+3,num2str(k),'color','b')
                        % rectangle('Position',[spots(k,1)-spots(k,6)/2,spots(k,2)-spots(k,7)/2,spots(k,6),spots(k,7)],'Curvature',[1,1],'EdgeColor','b')
                        %rectangle('Position',[spots(k,1)-spots(k,6),spots(k,2)-spots(k,7),spots(k,6)*2,spots(k,7)*2],'Curvature',[1,1],'EdgeColor','b')
                    end
                else
                for k=min(find(spots(:,9)==i)):max(find(spots(:,9)==i))
                    text(spots(k,1)+3,spots(k,2)+3,num2str(k),'color','b')
                    % rectangle('Position',[spots(k,1)-spots(k,6)/2,spots(k,2)-spots(k,7)/2,spots(k,6),spots(k,7)],'Curvature',[1,1],'EdgeColor','b')
                    rectangle('Position',[spots(k,1)-spots(k,6),spots(k,2)-spots(k,7),spots(k,6)*2,spots(k,7)*2],'Curvature',[1,1],'EdgeColor','b')
                end
                end
              
                hold off
                pause
            end
            
            
            %% LINK SPOTS INTO TRAJECTORIES
            if p.use_cursor==0
                %Start at 2nd frame
                if i>startFrame
                    %  if max(spots(:,9))==i
                    [spots]=LinkSpots4(spots, i, i-FrameInt, p);
                    %   end
                end
            end
            if p.show_text_output==1
                disp('trajectories determined')
            end
        else
            if p.show_text_output==1
                disp('No spots found')
            end
        end
    end
    %% FINAL PLOT
    if isempty(spots)==0
        if p.show_output==1
            if max(spots(:,10))>0
                subplot(2,3,2)
                imshow(image_data(:,:,startFrame),[],'InitialMagnification','fit')%HM modify magnification
                title('First frame with isolated spots superimposed')
                hold on
                plot(spots(spots(:,10)==0,1),spots(spots(:,10)==0,2),'o','color','r')
                subplot(2,3,1)
                imshow(image_data(:,:,startFrame),[],'InitialMagnification','fit')%HM modify magnification
                title('First frame with tracks superimposed')
                hold on
                for i=1:max(spots(:,10))
                    traj_color=rand(1,3);
                    subplot(2,3,1)
                    title('Trajectory intensities vs frame num')
                    if Ch==1
                        plot(spots(spots(:,10)==i,1),spots(spots(:,10)==i,2),'-o','color',traj_color)
                    elseif Ch==2
                        plot(spots(spots(:,10)==i,1)+round(size(image_data,2)/2+p.exclude_region),spots(spots(:,10)==i,2),'-o','color',traj_color)
                    end
                    subplot(2,3,3)
                    hold on
                    plot(spots(spots(:,10)==i,9),spots(spots(:,10)==i,5),'-o','color',traj_color)
                    
                end
                title('Trajectory intensities vs frame num')
                subplot(2,3,4)
                
                %    hist(spot_means)
                hist(spots(:,5))
                title('Histogram Trajectory Intensity Mean')
                subplot(2,3,5)
                
                %                 hist(spot_sd)
                %                 title('Histogram Trajectory Intensity SD')
            else
                if p.show_text_output==1
                    disp('No trajectories found')
                end
            end
        end
        %Assign spots to final variables
        if Ch==1
            if isempty(spots)==0
                SpotsCh1=spots;
            else
                SpotsCh1=0;
            end
        elseif Ch==2
            if isempty(spots)==0
                SpotsCh2=spots;
                % Transform Channel2 so spots are in the right place
                SpotsCh2(:,1)=SpotsCh2(:,1)+round(size(image_data,2)/2+p.exclude_region);
            else
                SpotsCh2=0;
            end
        end
        %  clear spots
    end
end
if p.show_text_output==0
 %   close(h)
end

end