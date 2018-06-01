p.noFrames=5; % number of frames to average over, starting from start_frame. (default: 5, end_frame-start_frame)

    p.subarray_halfwidth = 8; % (Default: 8 pixels). Halfwidth of image square subarray
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
    p.satPixelVal=10^10;
    
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
    p.show_output=1;
    % set this and graphs will appear for every spot!
    p.show_all_output=0;
    p.show_text_output=1;