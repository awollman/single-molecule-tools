%% initialise parameters
    params.pixelSize=0.08; %pixel size in microns
    params.frameLimitS=3; %num frames to use in stoich
    params.frameLimitD=3; %num frames to use in diffusion
    params.frameTime=0.005; % time between frames in seconds
    params.IsingleG=5000; %characteristic intensity of a single fluorophore
    params.IsingleR=5000; %characteristic intensity of a single fluorophore
    params.stoichMethod=1; %method for calculating stoichiomtry, see getStoichiometry
    params.bleachTime=5; %required for some stoich methods
    params.showOutput=1; %makes plots of each cell
    params.frameLimitAll=100; %number of frames in include in analysis
     params.transform=[-64,0]; %alignment between channels
    % maximum distance in pixels to link spots
    params.d=5;
    % maximum overlap integral to link spots
    params.overlap=0.75;
    % set if=0 link all spots regardless of frame
    %     if=1 link only spots in the same frame
    %     if=2 link spots in alternate frames (for ALEX)
    params.frameLinkMethod=0;
    
    % set so that spots are only assigned 1 partner
    params.Unique=1;
    
    %set to display graphs
    params.showOutput=1;

%% loop over data
sampleDir=uigetdir;
cd(sampleDir)
folderDir=dir('*cell*');
cellNo=1;
trackArrayCh1=[];
trackArrayCh2=[];

for c=1:length(folderDir)
    cd(folderDir(c).name)
    
    % look for files containging key words and load them in
    trackFile=dir('*TRACKS*');
    load(trackFile(1).name);
    segFile=dir('*segmentation*');
    load(segFile(1).name);
% segmentation might contain many cells so need to loop over each cell
    for s=1:size(CellObject,3)
    [trackArrayCh1Temp,trackArrayCh2Temp,SpotsCh1linked,SpotsCh2linked]=...
        colocalisedTrackAnalyser(SpotsCh1,SpotsCh2,CellObject(:,:,s),trackFile(1).name,cellNo,params);
    trackArrayCh1=cat(1,trackArrayCh1,trackArrayCh1Temp); %add all cells output into one array
    trackArrayCh2=cat(1,trackArrayCh2,trackArrayCh2Temp); %add all cells output into one array
    cellNo=cellNo+1;
    end
    cd ..
end
params.IgnoreEmpty=0; % will ignore empty cells if set to 1
h1=masterPlot(trackArrayCh1,params,2);
h2=masterPlot(trackArrayCh2,params,2);