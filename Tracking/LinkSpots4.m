function [spot_traj]=LinkSpots4(spots, framenum1, framenum0, p)
spot_traj=spots;
%d_01_max, Iratio_01_min, Iratio_01_max, SigmaRatio_01_min, SigmaRatio_01_max
% p.d_01_max, p.Iratio_01_min, p.Iratio_01_max, p.SigmaRatio_01_min, p.SigmaRatio_01_max
Spots2=spots(spots(:,9)==framenum1,:);
Spots1=spots(spots(:,9)==framenum0,:);

if isempty(Spots1) || isempty(Spots2)
else

 distanceMatrix=pdist2([Spots2(:,1),Spots2(:,2)],[Spots1(:,1),Spots1(:,2)]);
 %find all the best assignments for Spots2 in Spots1
 [SpotInd(:,1),SpotInd(:,2)]=min(distanceMatrix,[],2);
 SpotInd(:,3)=1:length(Spots2(:,1));
 % SpotInd contains distances, Spot1 indices, Spot2 indices
 % remove any too far away
 SpotInd(SpotInd(:,1)>p.d_01_max,:)=[];
  % remove any with unacceptable I or S ratios
SpotInd((p.Iratio_01_max > Spots1(SpotInd(:,2),5)./Spots2(SpotInd(:,3),5) > p.Iratio_01_min)==0,:)=[];
SpotInd((p.SigmaRatio_01_max > (Spots1(SpotInd(:,2),6)+Spots1(SpotInd(:,2),7))./(Spots2(SpotInd(:,3),6)+Spots2(SpotInd(:,3),7)) > p.SigmaRatio_01_min)==0,:)=[];
SpotIndSort=sortrows(SpotInd);
[~, ia, ~]=unique(SpotIndSort(:,2));
finalAssignments=SpotIndSort(ia,:);

% could  loop over non unique assignments in Spots1, then see if the orphan spot
% might link to something else

% give everything trajectory numbers
num4newTraj=length(Spots1(finalAssignments(Spots1(finalAssignments(:,2),10)<1,2),10));
if num4newTraj>0
Spots1(finalAssignments(Spots1(finalAssignments(:,2),10)<1,2),10)=(max(spots(:,10))+1):(max(spots(:,10))+num4newTraj);
end
Spots2(finalAssignments(:,3),10)=Spots1(finalAssignments(:,2),10);

% assign final variable
spot_traj(spot_traj(:,9)==framenum1,:)=Spots2;
spot_traj(spot_traj(:,9)==framenum0,:)=Spots1;
end
end