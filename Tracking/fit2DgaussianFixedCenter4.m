function [sdx, sdy, Icent, Rsquare] = fit2DgaussianFixedCenter4(image_data,Ibg_avg, Isp, x_centre, y_centre,myfit,options,p)

subarray_halfwidth=p.subarray_halfwidth;
    guess_sigma_Fit=p.guess_sigma_Fit;
    show_output=p.show_all_output;
sigmaFit_max=p.sigmaFit_max;
sigmaFit_min=p.sigmaFit_min;
if show_output==1
    pause on
end
image_data=double(image_data);
d = subarray_halfwidth;

% Initialise image array boundaries centred on spot centre
Iy_start=round(y_centre)-d;
Iy_end=round(y_centre)+d;
Ix_start=round(x_centre)-d;
Ix_end=round(x_centre)+d;

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

%Get the x and y co-ordinates of this grid
[Xs,Ys] = meshgrid(Ix_start:Ix_end, Iy_start:Iy_end);


% This is a matlab function that turns the x, y and z data within
% the matrix, into separate vectors:
[xData, yData, zData] = prepareSurfaceData(Xs, Ys, Idata);


% % Fit options:
options = fitoptions(myfit);
% 
% Starting parameter values:
% Icent_start = 2000;
sdx_start = guess_sigma_Fit;
sdy_start = guess_sigma_Fit;
% x0_start = x_centre;
% y0_start = y_centre;

%options.StartPoint = [Icent_start, sdx_start, sdy_start x0_start, y0_start];
options.StartPoint = [sdx_start, sdy_start];

% Lower parameter values:
% Icent_lower = 0;
sdx_lower = sigmaFit_min;
sdy_lower = sigmaFit_min;
% x0_lower = x_centre - disk_radius;
% y0_lower = y_centre - disk_radius;

%options.Lower = [Icent_lower, sdx_lower, sdy_lower x0_lower, y0_lower];
options.Lower = [sdx_lower, sdy_lower];

% Upper parameter values:plot(
% Icent_upper = inf;
sdx_upper = sigmaFit_max;
sdy_upper = sigmaFit_max;
% x0_upper = x_centre + disk_radius;
% y0_upper = y_centre + disk_radius;

%options.upper = [Icent_upper,  sdx_upper, sdy_upper, x0_upper, y0_upper];
options.upper = [sdx_upper, sdy_upper];






% Fitting of 2D gaussian to the data:
[fitobject, gof, output] = fit([xData, yData], zData, myfit, options, 'problem', {Ibg_avg,Isp, x_centre,y_centre});
if show_output==1
    fitdata=feval(fitobject,Xs,Ys);
    close all
    Isection1 = image_data(Iy_start:Iy_end,round(x_centre));
    subplot(1,3,1)
    plot(Isection1)
    hold on
 %   plot(fitdata(round(size(fitdata,1)*(x_centre/max(max(Xs)))),1:end),'r')
    plot(fitdata(1:end,round(size(fitdata,1)*((x_centre-min(min(Xs)))/(max(max(Xs))-min(min(Xs)))))),'r')
    hold off
    Isection2 = image_data(round(y_centre),Ix_start:Ix_end);
    subplot(1,3,2)
    plot(Isection2)
     hold on
  %  plot(fitdata(1:end,round(size(fitdata,1)*(y_centre/max(max(Ys))))),'r')
  %plot(fitdata(round(size(fitdata,1)*(y_centre/max(max(Ys)))),1:end),'r')
  plot(fitdata(round(size(fitdata,1)*((y_centre-min(min(Ys)))/(max(max(Ys))-min(min(Ys))))),1:end),'r')
    hold off
    subplot(1,3,3)
    plot(fitobject, [xData, yData], zData);
    pause
end
output_params = coeffvalues(fitobject);
sdx=output_params(1);
sdy=output_params(2);
Icent=feval(fitobject, [x_centre, y_centre]);
Rsquare=gof.rsquare;
end