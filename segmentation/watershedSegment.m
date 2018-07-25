function NewMask=watershedSegment(ImageFrame,Mask,seedMask)
oldMask=Mask;
I_eq_c = imcomplement(ImageFrame); %complement the image because we are about to apply the watershed transform, which identifies low points, not high points.
mask_em=seedMask;
mask_em(mask_em>1)=1;

mask_em(1,:)=0;
mask_em(:,1)=0;
mask_em(end,:)=0;
mask_em(:,end)=0;

%figure; imshow(mask_em,[])
Mask=double(Mask)+double(seedMask);
Mask(Mask>1)=1;
Mask(1,:)=0;
Mask(:,1)=0;
Mask(end,:)=0;
Mask(:,end)=0;
%figure; imshow(Mask,[])
I_mod = imimposemin(I_eq_c, ~Mask | mask_em); %modify image so background and extended maxima pixels forced to be only local minima in image.
L_cells = watershed(I_mod);
%NewMask=L_cells;
% only write the new mask if watershedding has worked

% rgb_cells = label2rgb(L_cells,'jet',[.5 .5 .5]);
% figure;
% imshow(rgb_cells,'InitialMagnification','fit')
% title('Watershed transform of D')
%NewMask=L_cells;
if range(range(L_cells))>0
    NewMaskTemp=zeros(size(ImageFrame,1),size(ImageFrame,2),max(L_cells(:))-1);
    
    for lbl = 1:max(L_cells(:))
        NewMaskTemp(:,:,lbl) = L_cells==lbl;
    end
    NewMask=NewMaskTemp;
else
    NewMask=oldMask;
end



% for lbl = 1:max(L_cells(:))
%     OverlapMask=seedMask(:,:,l).*NewMaskTemp(:,:,lbl);
%     if sum(OverlapMask(:))>0
%         NewMask(:,:,lbl)=NewMaskTemp(:,:,lbl);
%     end
% end


end