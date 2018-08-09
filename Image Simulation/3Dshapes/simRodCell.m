function projection_image=simRodCellShell(radius, lengthVar, realCenter, orientation,psf,sizeImage)

Center=[sizeImage(2)/2,sizeImage(1)/2];
projection_imageCylinder=simCylinder(radius, lengthVar, Center, psf,sizeImage);
Center1=Center;
Center1(2)=Center(2)-lengthVar;
projection_imageSphere1=simHalfSphere(radius, Center, psf, size(projection_imageCylinder));
sphereShift=circshift(projection_imageSphere1,[-lengthVar/2,0]);

%projection_imageSphere1=zeros(size(projection_imageCylinder));
projection_imageTemp=projection_imageCylinder+sphereShift+circshift(flipud(sphereShift),[1,0]);
projection_imageRot=imrotate(projection_imageTemp,orientation,'nearest','crop');
projection_image=circshift(projection_imageRot,[realCenter(1)-Center(1),realCenter(2)-Center(2)]);
end