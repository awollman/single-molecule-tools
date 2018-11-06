% Outputs:
% NucSpotTotCorr=total number of molecules in spots in the nucleus
% CellSpotTotCorr=total number of molecules in spots in the cytoplasm
% N_nuc_poolCorr=total number of molecules in the pool in the nucleus
% N_cyto_poolCorr=total number of molecules in the pool in the cytoplasm
% Cell_Rad=cell radius
% Nuc_Rad=nucleus radius
% Inputs
% ImageData=single fluorescence image frame
% AutoFluorData=either an autofluorescent image or the mean autofluorescent
% background
% CellSeg=cell segmentation 1 CELL ONLY
% NucSeg=nucleus segmentation
% Spots=Spots array from ADEMScode only for spots in the frame
% Isingle=characteristic intensity of fluorophore
% psf=model psf for deconvolution, can use the one in theoreticalPSF.mat

function [NucSpotTotCorr, CellSpotTotCorr, N_nuc_poolCorr, N_cyto_poolCorr, Cell_Rad, Nuc_Rad]=...
    CopyNumberSphereCell(imageData,AutoFluorData, CellSeg, NucSeg, Spots, Isingle,psf,params)

if nargin<8
    params.spotAvoid=0;
    params.showOutput=1;
    params.inner_circle_radius=5;
    params.method=1;
    params.depthfield=3;
end

%% Generate spot matrix to avoid areas where there are spots
se=strel('disk',params.inner_circle_radius);
SpotMatrix=zeros(size(NucSeg));
for i=1:size(Spots,1)
    SpotMatrix(round(Spots(i,2)),round(Spots(i,1)))=1;
end
SpotMatrixDisk=imdilate(SpotMatrix,se);

spotAvoid=ones(size(NucSeg))-SpotMatrixDisk; 
%% Cytoplasmic concentration
se2=strel('disk',2);
CellPixels=(CellSeg-NucSeg);
CellPixels(CellPixels<0)=0;
CellAvoidPixels=CellSeg-imdilate(NucSeg,se);
CellAvoidPixels(CellAvoidPixels<0)=0;
CellNoSpots=imerode(CellSeg,se2).*spotAvoid;
J=AutoFluorData;
BGcorrImage=double(imageData)-double(J);
Cell_stats=regionprops(CellSeg,'Area','Centroid','MinorAxisLength');

Cell_area=Cell_stats.Area;
Cell_centroid=Cell_stats.Centroid;
Cell_Rad=(Cell_area/pi).^0.5;
%     Cell_Rad(cell_counter)=Cell_stats.MinorAxisLength/2;
projection_image=simSphere(Cell_Rad, Cell_centroid, psf,size(imageData));
%    imshow(BGcorrImage,[])
Nuc_stats=regionprops(NucSeg,'Area','Centroid');
Nuc_centroid=Nuc_stats.Centroid;
Nuc_area=Nuc_stats.Area;
Nuc_Rad=(Nuc_area/pi).^0.5;
Nuc_projection_image=simSphere(Nuc_Rad, Nuc_centroid, psf,size(imageData));
cyto_projection_image=projection_image-Nuc_projection_image;

Rp=BGcorrImage(CellNoSpots==1)./Isingle;
Cc=cyto_projection_image(CellNoSpots==1);
Cn=Nuc_projection_image(CellNoSpots==1);

A=[Cc,Cn];
lb=[0,0];
ub=[Inf,Inf];
X=lsqlin(A,Rp,[],[],[],[],lb,ub,[]);
CytoConc= X(1);
NucConc= X(2);



NucSpotTot=0;
CellSpotTot=0;
for l=1:size(Spots,1)
    
    if NucSeg(round(Spots(l,2)),round(Spots(l,1)))>0
        
        NucSpotTot=NucSpotTot+Spots(l,5)/Isingle;
        
    else
        if CellSeg(round(Spots(l,2)),round(Spots(l,1)))>0
            CellSpotTot=CellSpotTot+Spots(l,5)/Isingle;
            
        end
        
    end
    
end

CellVol=4/3.*pi.*Cell_Rad.^3;
NucVol=4/3.*pi.*Nuc_Rad.^3;
%NucVol=4/3.*pi.*((Nuc_area)./pi).^1.5;
if CellVol>NucVol
CytoVol=CellVol-NucVol;
else
    CellVol=4/3.*pi.*30.^3;
NucVol=4/3.*pi.*13.^3;
end

N_cyto_poolTemp=CytoConc.*CytoVol;
N_nuc_poolTemp=NucConc*NucVol;

if  params.method==1
    
    N_cyto_pool=N_cyto_poolTemp;
N_nuc_pool=N_nuc_poolTemp;
    
else % pin the total copy number to the summed intensity corrected for diffraction
totalProj=cyto_projection_image+Nuc_projection_image;
%sum(sum(totalProj(:))/sum(totalProj(CellSeg==1)))
totalCopyNumber=mean(Rp)*sum(CellSeg(:))*sum(sum(totalProj(:))/sum(totalProj(CellSeg==1)));

 N_cyto_pool=totalCopyNumber*(N_cyto_poolTemp/(N_cyto_poolTemp+N_nuc_poolTemp));
  N_nuc_pool=totalCopyNumber*(N_nuc_poolTemp/(N_cyto_poolTemp+N_nuc_poolTemp));
end

NucSpotTotCorr=(NucSpotTot.*NucVol./(Nuc_area*2*params.depthfield));
CellSpotTotCorr=(CellSpotTot.*CytoVol./((Cell_area-Nuc_area)*2*params.depthfield));
N_nuc_poolCorr=(N_nuc_pool-(NucSpotTotCorr-NucSpotTot));
N_cyto_poolCorr=(N_cyto_pool-(CellSpotTotCorr-CellSpotTot));


if params.showOutput==1
    figure;
    subplot(1,3,1)
    imshow(imageData,[])
    hold on
    [row,col]=find(bwperim(CellSeg));
    scatter(col,row,20,'filled')
        [row,col]=find(bwperim(NucSeg));
    scatter(col,row,20,'filled')
    title('Original Data and segmentation')
    subplot(1,3,2)
    imshow(double(imageData).*double(spotAvoid),[])
    title('Intensity map used for copy number determination')
    %colorbar
    subplot(1,3,3)
    imshow(cat(3,mat2gray(Nuc_projection_image),mat2gray(cyto_projection_image),zeros(size(cyto_projection_image))))
        title('Simulated image used to deconvolve, Nuc=red, cyto=green')

end


end