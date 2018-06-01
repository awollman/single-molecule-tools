function [sdx, sdy, Icent] = iterate1DgaussianFixedCenter3(image_data,Ibg_avg, Isp, x_centre, y_centre,p)

subarray_halfwidth=p.subarray_halfwidth;
    guess_sigma_Fit=p.guess_sigma_Fit;
    show_output=p.show_all_output;
image_data=double(image_data);

d = subarray_halfwidth;

% Initialise image array boundaries centred on spot centre
Iy_start=round(y_centre)-d;
Iy_end=round(y_centre)+d;
Ix_start=round(x_centre)-d;
Ix_end=round(x_centre)+d;

% NewYcent=y_centre-Iy_start
% 
% NewXcent=x_centre-Ix_start

% If the boundaries are outside image boundary, move them
if (round(y_centre)-d)<1
    
    Iy_start=1;
    Iy_end=2*d;
end
if (round(y_centre)+d)>size(image_data,1)
    Iy_start=size(image_data,1)-2*d;
    Iy_end=size(image_data,1);
end
if (round(x_centre)-d)<1
    Ix_start=1;
    Ix_end=2*d;
end
if (round(x_centre)+d)>size(image_data,2)
    Ix_start=size(image_data,2)-2*d;
    Ix_end=size(image_data,2);
end

%Define a grid around found spot centres to fit gaussian to.
Idata = image_data(Iy_start:Iy_end,Ix_start:Ix_end);

NewYcent=y_centre-Iy_start;
% 
NewXcent=x_centre-Ix_start;


%Get the x and y co-ordinates of this grid
%[Xpos,Ypos] = meshgrid(Ix_start:Ix_end, Iy_start:Iy_end);
Xsize=d*2+1;
Ysize=d*2+1;
[Xpos,Ypos] = meshgrid(1:size(Idata,2),1:size(Idata,1));
error_set=0.005;

PSFwidthInitial=guess_sigma_Fit/2;
k=1;
%Set intial residuals to be massive
Residuals(k)=10^12;
changeInt=0.99;
switchNo=1;
PSFwidth2=PSFwidthInitial;
while (changeInt>error_set && PSFwidth2>p.sigmaFit_min && PSFwidth2<p.sigmaFit_max)
 %   while (k<10)
if k==1
    PSFwidth=PSFwidthInitial;
    else
    PSFwidth=PSFwidth*(1+switchNo*changeInt);
end
PSFwidth2=PSFwidth;
%GaussFrame=Isp./(2.*pi.*PSFwidth.*PSFwidth)*exp(-(((Xpos-size(Idata,2)/2).^2)./(2.*PSFwidth^2)+((Ypos-size(Idata,1)/2).^2)./(2.*PSFwidth^2)));
% subplot(2,2,1)
% surf(GaussFrame)
% subplot(2,2,2)
% surf(Idata)
% subplot(2,2,3)
% surf(Xpos)
% subplot(2,2,4)
% surf(Ypos)
GaussFrame=Isp./(2.*pi.*PSFwidth.*PSFwidth)*exp(-(((Xpos-NewYcent).^2)./(2.*PSFwidth^2)+((Ypos-NewXcent).^2)./(2.*PSFwidth^2)));


ResidualImage=(Idata-Ibg_avg-GaussFrame);
Residuals(k+1)=sum(abs(ResidualImage(1:end)));

if Residuals(k+1)<Residuals(k)
   changeInt=changeInt/1.5;
else
    switchNo=-switchNo;
end

%figure; surf(Idata-Ibg_avg-GaussFrame)

k=k+1;
if k>20
    break
end
end
if PSFwidth<p.sigmaFit_max
sdx=PSFwidth;
sdy=PSFwidth;
else
    sdx=p.sigmaFit_max;
sdy=p.sigmaFit_max;
end
Icent=GaussFrame(d,d);

if show_output==1
subplot(2,2,1)
plot(Residuals(2:end))
xlabel('iterations')
ylabel('Residuals')
subplot(2,2,2)
[xData, yData, zData] = prepareSurfaceData(Xpos, Ypos, GaussFrame);
surf(Idata-Ibg_avg)
hold on
scatter3(xData, yData, zData,100,'r','.')
zlim([0 max(zData)])
subplot(2,2,3)
plot(Idata(:,8)-Ibg_avg,'.')
hold on
plot(GaussFrame(:,8))
subplot(2,2,4)
plot(Idata(8,:)-Ibg_avg,'.')
hold on
plot(GaussFrame(8,:))
pause
close
end