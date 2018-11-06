


%% loop over data
sampleDir=uigetdir;
cd(sampleDir)
folderDir=dir('*cell*');
cellNo=1;
allSpots=[];
for c=1:length(folderDir)
    cd(folderDir(c).name)

    % look for files containging key words and load them in
    trackFile=dir('*TRACKS*');
    load(trackFile(1).name);
    allSpots=cat(1,allSpots,SpotsCh1);
    cd ..
end

figure;
[counts,x]=ksdensity(allSpots(:,5));
plot(x,counts)
xlabel('Intensity')
ylabel('Probability')
