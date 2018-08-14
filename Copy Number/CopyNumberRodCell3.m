function [poolNo, CellSpotTot, noSpots, concentration, area, volume]=CopyNumberRodCell3(imageData, segmentation, BGav, SpotsCh1,psf, Isingle,params)
if nargin<7
    params.spotAvoid=0;
    params.showOutput=1;
    params.inner_circle_radius=5;
    params.method=1;
end

 [cellCoord(:,2), cellCoord(:,1)]=find(segmentation);
spotInd=ismember(round(SpotsCh1(:,1:2)),cellCoord,'rows');
Spots=SpotsCh1(spotInd,:);
%% avoid foci within the depth of field
if params.spotAvoid==1
    se=strel('disk', params.inner_circle_radius);
    SpotMatrix=zeros(size(imageData));
    for k=1:size(Spots,1)
        SpotMatrix(round(Spots(k,2)),round(Spots(k,1)))=1;
    end
    SpotMatrixDisk=imdilate(SpotMatrix,se);
    SpotAvoid=ones(size(imageData))-SpotMatrixDisk;
    CellIntensities=double(imageData).*segmentation.*SpotAvoid;
else
    CellIntensities=double(imageData).*segmentation;
end
% I_cell_median=median(CellIntensities(CellIntensities>0))-BGav;
% I_cell_std=std(CellIntensities(CellIntensities>0));

stats=regionprops(segmentation,'Centroid','MajorAxisLength','MinorAxisLength','Orientation','Area');
area=stats.Area;
Centroid=[stats.Centroid(1),stats.Centroid(2)];
projection_image=simRodCell(round(stats.MinorAxisLength/2),round(stats.MajorAxisLength),round(Centroid) ,...
    90+stats.Orientation, psf, size(imageData));

CellIntensities(CellIntensities>0)=CellIntensities(CellIntensities>0)-BGav;
volume = (pi*(stats.MinorAxisLength/2)^2 * (stats.MajorAxisLength-(stats.MinorAxisLength/2))) + (4*pi/3 * (stats.MinorAxisLength/2)^3);

switch params.method
    case 1 % solve linear equations
        concentration=(projection_image(CellIntensities~=0).*Isingle)\CellIntensities(CellIntensities~=0);
        poolNo=volume*concentration;

    case 2 % sum intensity in image and correct for PSF convolution
       correctionFactor=sum(projection_image)/sum(projection_image.*segmentation)
        poolNo=(mean(CellIntensities(CellIntensities~=0))/Isingle)*area*correctionFactor;
        concentration=poolNo/volume;
end


% calculate num in spots
CellSpotTot=0;
noSpots=0;
for l=1:size(Spots,1)
    
        if segmentation(round(Spots(l,2)),round(Spots(l,1)))>0
            CellSpotTot=CellSpotTot+Spots(l,5)/Isingle;
            noSpots=noSpots+1;
            
        end
    
end

%% add in spot correction

if params.showOutput==1
    figure;
    subplot(2,2,1)
    imshow(imageData,[])
    hold on
    [row,col]=find(bwperim(segmentation));
    scatter(col,row,20,'filled')
    title('Original Data and segmentation')
    subplot(2,2,2)
    imshow(CellIntensities,[])
    title('Intensity map used for copy number determination')
    colorbar
    subplot(2,2,3)
    imshow(projection_image,[])
        title('Simulated image used to deconvolve')
    colorbar
        subplot(2,2,4)
    imshow(CellIntensities/(Isingle*concentration),[])
            title('Concentration Map')
    colorbar
end

end