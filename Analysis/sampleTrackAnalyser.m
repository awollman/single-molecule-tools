%% initialise parameters
    params.pixelSize=0.08; %pixel size in microns
    params.frameLimitS=4; %num frames to use in stoich
    params.frameLimitD=4; %num frames to use in diffusion
    params.frameTime=0.005; % time between frames in seconds
    params.Isingle=5000; %characteristic intensity of a single fluorophore
    params.stoichMethod=3; %method for calculating stoichiomtry, see getStoichiometry
    params.bleachTime=5; %required for some stoich methods
    params.showOutput=1; %makes plots of each cell
    params.frameLimitAll=10; %number of frames in include in analysis



%% loop over data
sampleDir=uigetdir;
cd(sampleDir)
folderDir=dir('*cell*');
cellNo=1;
trackArray=[];
allSpots=[];
for c=1:length(folderDir)
    cd(folderDir(c).name)
    
    % look for files containging key words and load them in
    trackFile=dir('*TRACKS*');
    load(trackFile(1).name);
    segFile=dir('*segmentation*');
    load(segFile(1).name);
    if params.showOutput==1
        h1=figure;
        imshow(frame_average,[])
        hold on
    end
% segmentation might contain many cells so need to loop over each cell
    for s=1:size(CellObject,3)
    [trackArrayTemp,spotsInTracks]=trackAnalyser(SpotsCh1,CellObject(:,:,s),trackFile(1).name,cellNo,params);
    spotsInTracks(:,13)=cellNo;
    allSpots=cat(1,allSpots,spotsInTracks);
    trackArray=cat(1,trackArray,trackArrayTemp); %add all cells output into one array
    if params.showOutput==1
        figure(h1)
        [row,col]=find(bwperim(CellObject(:,:,s)));
        rp=regionprops(CellObject(:,:,s),'centroid');
        scatter(col, row)
        text(rp.Centroid(1), rp.Centroid(2),num2str(cellNo),'color','red')
    end
    cellNo=cellNo+1;
    end
    cd ..
end
params.IgnoreEmpty=0; % will ignore empty cells if set to 1
h1=masterPlot(trackArray,params,1);