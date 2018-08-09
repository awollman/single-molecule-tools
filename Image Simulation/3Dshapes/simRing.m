% creates the projection of a ring at Center, radius and width rotated
% 90 degrees in 3 permutations

function projection_image=PSFintRing3(radius, Center, width, permutation, psf, sizeImage)

%width=3;
image3D=zeros(sizeImage,sizeImage,sizeImage);
%Center=[(sizeImage+1)/2,(sizeImage+1)/2,(sizeImage+1)/2];
zwidth=width*0.5;
switch permutation
    case 1
        for i=1:sizeImage
            for j=1:sizeImage
                for k=1:sizeImage
                    radialDist=(Center(1)-i)^2+(Center(2)-j)^2+(Center(3)-k)^2;
                    if abs((Center(3)-k))<zwidth
                        if radialDist<((radius)^2) && radialDist>((radius-width)^2)
                            image3D(i,j,k)=1;
                        end
                    end
                end
            end
        end
        image3Drotate=image3D;
    case 2
        for i=1:sizeImage
            for j=1:sizeImage
                for k=1:sizeImage
                    radialDist=(Center(1)-i)^2+(Center(2)-j)^2+(Center(3)-k)^2;
                    if abs((Center(2)-j))<zwidth
                        if radialDist<((radius)^2) && radialDist>((radius-width)^2)
                            image3D(i,j,k)=1;
                        end
                    end
                end
            end
        end
        image3Drotate=image3D;
    case 3
        for i=1:sizeImage
            for j=1:sizeImage
                for k=1:sizeImage
                    radialDist=(Center(1)-i)^2+(Center(2)-j)^2+(Center(3)-k)^2;
                    if abs((Center(1)-i))<zwidth
                        if radialDist<((radius)^2) && radialDist>((radius-width)^2)
                            image3D(i,j,k)=1;
                        end
                    end
                end
            end
        end
        image3Drotate=image3D;
end

image=zeros(sizeImage,sizeImage);
projection_image=zeros(sizeImage,sizeImage);
for i=1:size(image,2)
    for j=1:size(image,1)
        for k=1:sizeImage
            if image3Drotate(i,j,k)==1
                
                %     projection_image(j,i)=projection_image(j,i)+1;
                projection_image=projection_image+...
                    psf((128-j+1):(128+size(image,1)-j),...
                    (128-i+1):(128+size(image,2)-i),k);
                %              mesh( psf((128-i+1):(256-i),...
                %                     (160-j+1):(224-j),k))
                %                 pause
                
            end
        end
    end
end

%mesh(projection_image)
%imshow(projection_image,[])
end