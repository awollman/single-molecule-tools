function [xNew,yNew]=MergeCoincidentCandidates2(x_estimate,y_estimate, d_min)
i=1;
spots_temp=zeros(size(x_estimate,1),12);
spots_temp(:,1)=x_estimate;
spots_temp(:,2)=y_estimate;
spots_temp(:,9)=i;
[spots2, new_spots]=MergeCoincidentSpots3(spots_temp, i, d_min);

 spots2(spots2(:,1)==100000,:)=[];
xNew=spots2(:,1);
yNew=spots2(:,2);
end