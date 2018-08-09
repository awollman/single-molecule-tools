function projection_image=simCylinderShell(radius, lengthVar, width,Center, psf,sizeImage)
%radius=30;
%Center=[32, 64];
%psf;
maxZ=max([sizeImage(1),sizeImage(2),radius]);
projection_image=zeros(sizeImage);
for i=1:sizeImage(2)
    for j=1:sizeImage(1)
        for k=1:maxZ
            if (Center(1)-i)^2+(33-k)^2<radius^2
                if (Center(1)-i)^2+(33-k)^2>((radius-width)^2)
                    if abs((Center(2)-j)) < lengthVar/2
                        %  projection_image(j,i)=projection_image(j,i)+1;
               projection_image=projection_image+...
            psf((size(psf,1)/2-j+1):(size(psf,1)/2+sizeImage(1)-j),...
                   (size(psf,2)/2-i+1):(size(psf,2)/2+sizeImage(2)-i),k);
                        %              mesh( psf((128-i+1):(256-i),...
                        %                     (160-j+1):(224-j),k))
                        %                 pause
                    end
                end
            end
        end
    end
end

%mesh(projection_image)
end