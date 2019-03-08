


%% loop over data
sampleDir=uigetdir;
cd(sampleDir)
folderDir=dir('*e*');
cellNo=1;
allSpots=[];
for c=1:length(folderDir)
    cd(folderDir(c).name)
    cellCoord=[];
%try
    % look for files containging key words and load them in
    trackFile=dir('*TRACKS*');
    load(trackFile(1).name);
        segFile=dir('*segmentation*');
    load(segFile(1).name);
    
    [cellCoord(:,2), cellCoord(:,1)]=find(sum(Mask,3));
    SpotsCh2(:,1)=SpotsCh2(:,1)-500;
spotInd=ismember(round(SpotsCh2(:,1:2)),cellCoord,'rows');
spots=SpotsCh2(spotInd,:);
    allSpots=cat(1,allSpots,spots);
% catch
% end
    cd ..
end

figure;
%[counts,x]=ksdensity(allSpots(:,5),'bandwidth',1);
[counts,x]=ksdensity(allSpots(:,5));

plot(x,counts)
xlabel('Intensity')
ylabel('Probability')
