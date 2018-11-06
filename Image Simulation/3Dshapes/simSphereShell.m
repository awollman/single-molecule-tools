function projection_image=simSphereShell(radius, Center, width, psf, sizeImage)
maxZ=size(psf,3);
projection_image=zeros(sizeImage);
if size(psf,1)<2*sizeImage(1)
    psf= padarray(psf,[round((2*sizeImage(1)-size(psf,1))/2),0],'both');
end

if size(psf,2)<2*sizeImage(2)
    psf= padarray(psf,[0,round((2*sizeImage(2)-size(psf,2))/2)],'both');
end

for i=1:sizeImage(2)
    for j=1:sizeImage(1)
        for k=1:maxZ
            if (Center(1)-i)^2+(Center(2)-j)^2+(ceil(maxZ/2)-k)^2<radius^2
                if (Center(1)-i)^2+(Center(2)-j)^2+(ceil(maxZ/2)-k)^2>((radius-width)^2)
           %     projection_image(j,i)=projection_image(j,i)+1;

               projection_image=projection_image+...
            psf((size(psf,1)/2-j+1):(size(psf,1)/2+sizeImage(1)-j),...
                   (size(psf,2)/2-i+1):(size(psf,2)/2+sizeImage(2)-i),k);

                end
            end
        end
    end
end
end
