function [trackArrayCh1,trackArrayCh2,SpotsCh1linked,SpotsCh2linked]=colocalisedTrackAnalyser(SpotsCh1,SpotsCh2,segmentedCell,tracksFile,cellNo,params)    
if nargin<6
    params.pixelSize=0.08; %pixel size in microns
    params.frameLimitS=4; %num frames to use in stoich
    params.frameLimitD=4; %num frames to use in diffusion
    params.frameTime=0.005; % time between frames in seconds
    params.IsingleG=5000; %characteristic intensity of a single fluorophore
        params.IsingleR=5000; %characteristic intensity of a single fluorophore
    params.stoichMethod=1; %method for calculating stoichiomtry, see getStoichiometry
    params.bleachTime=5; %required for some stoich methods
    params.showOutput=1; %makes plots of each cell
    params.frameLimitAll=10; %number of frames in include in analysis
     params.transform=[0,0]; %alignment between channels
    % maximum distance in pixels to link spots
    params.d=5;
    % maximum overlap integral to link spots
    params.overlap=0.75;
    
    % set if=0 link all spots regardless of frame
    %     if=1 link only spots in the same frame
    %     if=2 link spots in alternate frames (for ALEX)
    params.frameLinkMethod=1;
    
    % set so that spots are only assigned 1 partner
    params.Unique=1;
    
    %set to display graphs
    params.showOutput=1;
end

SpotsCh2(:,1)=SpotsCh2(:,1)+params.transform(1);
SpotsCh2(:,2)=SpotsCh2(:,2)+params.transform(2);
params.Isingle=params.IsingleG;
[trackArrayCh1,spotsInTracksCh1]=trackAnalyser(SpotsCh1,segmentedCell,tracksFile,cellNo,params);
params.Isingle=params.IsingleR;
[trackArrayCh2,spotsInTracksCh2]=trackAnalyser(SpotsCh2,segmentedCell,tracksFile,cellNo,params);
params.transform=[0,0]; %need to reset this so Colocaliser2 does not also transform coords
if ~isempty(spotsInTracksCh1) && ~isempty(spotsInTracksCh2)
[SpotsCh1linked, SpotsCh2linked]=Colocaliser2(spotsInTracksCh1,spotsInTracksCh2,params);

    
 % this is section is a loop for now, possible with arrays but
            % difficult, safer to stick to loop for now
            [ch1Traj,ch1Ind]=unique(SpotsCh1linked(SpotsCh1linked(:,14)>0,10)); % linked trajectories in ch1
            ch2Traj=SpotsCh1linked(SpotsCh1linked(:,14)>0,15);
            for tr=1:length(ch1Traj) %loop over linked trajectory numbers
                % link tractory in ch1 with the longest linked traj in ch2
                ch2Traj(tr)=mode(SpotsCh1linked(SpotsCh1linked(:,14)>0 & SpotsCh1linked(:,10)==ch1Traj(tr),15));
                trackArrayCh1(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr),8)=ch2Traj(tr);
                % length of time in frames linked to a trajectory
                trackArrayCh1(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr),9)=sum(SpotsCh1linked(:,10)==ch1Traj(tr) & SpotsCh1linked(:,15)==ch2Traj(tr));
                
                trackArrayCh1(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr),10)=find(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr));
                % distance from other traj
                trackArrayCh1(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr),11)=mean(SpotsCh1linked(SpotsCh1linked(:,10)==ch1Traj(tr) & SpotsCh1linked(:,15)==ch2Traj(tr),16));
                % link corresponding tractory in ch2 with ch1
                trackArrayCh2(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr),8)=ch1Traj(tr);
                % assigne the same link time to ch2
                trackArrayCh2(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr),9)=sum(SpotsCh1linked(:,10)==ch1Traj(tr) & SpotsCh1linked(:,15)==ch2Traj(tr));
                trackArrayCh2(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr),10)=find(trackArrayCh1(:,7)==trackArrayCh1(end,7) & trackArrayCh1(:,6)==ch1Traj(tr));
                
                trackArrayCh2(trackArrayCh2(:,7)==trackArrayCh2(end,7) & trackArrayCh2(:,6)==ch2Traj(tr),11)=mean(SpotsCh1linked(SpotsCh1linked(:,10)==ch1Traj(tr) & SpotsCh1linked(:,15)==ch2Traj(tr),16));
            end
else
    SpotsCh1linked=SpotsCh1;
    SpotsCh2linked=SpotsCh2;
    
end
end