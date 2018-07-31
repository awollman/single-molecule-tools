function [trackArrayCh1,trackArrayCh2,SpotsCh1linked,SpotsCh2linked]=colocalisedTrackAnalyser(SpotsCh1,SpotsCh2,segmentedCell,tracksFile,cellNo,params)    
[trackArrayCh1,spotsInTracksCh1]=trackAnalyser(SpotsCh1,segmentedCell,tracksFile,cellNo,params);
[trackArrayCh2,spotsInTracksCh2]=trackAnalyser(SpotsCh2,segmentedCell,tracksFile,cellNo,params);

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
end