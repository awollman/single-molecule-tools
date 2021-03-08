function [spots_temp2,Spot1indN]=MergeCoincidentSpots4(spots_temp, d_min)
spots_temp2=spots_temp;
 distanceList=pdist([spots_temp(:,1),spots_temp(:,2)]);
 distanceMatrix=squareform(distanceList);
 removeList=tril((distanceMatrix<d_min)-eye(length(spots_temp(:,1)))); %need to remove top half of the triangle
  [Spot1ind,Spot2ind] = find(removeList);
  [Spot1indN,Spot2indN] = find(~removeList);
  for sp=1:length(Spot1ind)
      spots_temp2(Spot1ind(sp),1)=mean([spots_temp(Spot1ind(sp),1),spots_temp(Spot2ind(sp),1)]);
      spots_temp2(Spot1ind(sp),2)=mean([spots_temp(Spot1ind(sp),2),spots_temp(Spot2ind(sp),2)]);
  end
  spots_temp2(Spot2ind,:)=[];
% 
 
end

%distanceMatrix=pdist([SpotsCh1(:,1),SpotsCh1(:,2)])