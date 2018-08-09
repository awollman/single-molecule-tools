function projection_image=simHalfSphere(radius, Center, psf, sizeImage)
%radius=30;
%Center=[32, 64];
%psf;
maxZ=max([sizeImage(1),sizeImage(2),radius]);
projection_image=zeros(sizeImage);
for i=1:sizeImage(2)
    for j=1:sizeImage(1)
        for k=1:maxZ
            if (Center(1)-i)^2+(Center(2)-j)^2+(ceil(maxZ/2)-k)^2<radius^2
                %  if (Center(1)-i)^2+(Center(2)-j)^2+(33-k)^2>((radius-width)^2)
                if (Center(2)-j)>0
                    %     projection_image(j,i)=projection_image(j,i)+1;
                    
                    projection_image=projection_image+...
                        psf((size(psf,1)/2-j+1):(size(psf,1)/2+sizeImage(1)-j),...
                        (size(psf,2)/2-i+1):(size(psf,2)/2+sizeImage(2)-i),k);
                end
                %   end
            end
        end
    end
end

%mesh(projection_image)
end