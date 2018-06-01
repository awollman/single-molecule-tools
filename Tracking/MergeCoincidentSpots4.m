function [spots_temp2]=MergeCoincidentSpots4(spots_temp, d_min)
 distanceList=pdist([spots_temp(:,1),spots_temp(:,2)]);
 distanceMatrix=squareform(distanceList);
 removeList=(distanceMatrix<d_min)-eye(length(spots_temp(:,1)));
  [Spot1ind,Spot2ind] = find(removeList);
% 
 
end

%distanceMatrix=pdist([SpotsCh1(:,1),SpotsCh1(:,2)])