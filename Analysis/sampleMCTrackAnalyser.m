%% Multi-compartment track analyser

%Example script to loop over tracked and segmented data and analyse
%trajectories. Also categorises tracks into multiple compartments.

% Requires TRACKS file with frame_average image and SpotsCh1
% SEGMENTATION for the cell called CellObject, N x M x C = M,N pixels, C
% cells
% COMPARTMENTSEGMENTATION, same as for the cells


%% initialise parameters
    params.pixelSize=0.05; %pixel size in microns
    params.frameLimitS=4; %num frames to use in stoich
    params.frameLimitD=4; %num frames to use in diffusion
    params.frameTime=0.005; % time between frames in seconds
    params.Isingle=191; %characteristic intensity of a single fluorophore
    params.stoichMethod=3; %method for calculating stoichiomtry, see getStoichiometry
    params.bleachTime=5; %required for some stoich methods
    params.showOutput=1; %makes plots of each cell
    params.frameLimitAll=20; %number of frames to include in analysis



%% loop over data
sampleDir=uigetdir;
cd(sampleDir)
folderDir=dir('*Video*');
cellNo=1;
trackArray=[];
for c=1:length(folderDir)
    cd(folderDir(c).name)
    
    % look for files containging key words and load them in
    
    trackFile=dir('*TRACK*');
    load(trackFile(1).name);
        segFile=dir('*segmentation*');
    load(segFile(1).name);
    segFile2=dir('*compartmentSegmentation*');
    load(segFile2(1).name);
% segmentation might contain many cells so need to loop over each cell
    for s=1:size(CellObject,3)
        finalMask(:,:,1)=CellObject(:,:,s);
        for a=1:size(compMask)
            if sum(sum(CellObject(:,:,s).*compMask(:,:,a)))>0
                finalMask(:,:,end+1)=compMask(:,:,a);
            end
        end
    [trackArrayTemp,spotsInTracks]=trackAnalyser(SpotsCh1,finalMask,trackFile(1).name,cellNo,params);
    trackArray=cat(1,trackArray,trackArrayTemp); %add all cells output into one array
    cellNo=cellNo+1;
    end
    cd ..
end
params.IgnoreEmpty=0; % will ignore empty cells if set to 1
h1=masterPlot(trackArray,params,1);
