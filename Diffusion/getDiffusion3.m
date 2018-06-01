function [DiffusionConst,MSD,tau,LocPrecision]=getDiffusion3(Spot,frameTime,pixelSize,nPoints)

       G = [497.4  ,-54.99, 1050];
       b = [ 5.056  ,-1.184, 11.3];
       x=[mean(Spot(1:nPoints,11)),mean(Spot(1:nPoints,11))-std(Spot(1:nPoints,11)),mean(Spot(1:nPoints,11))+std(Spot(1:nPoints,11))];
       if x(2)<0
           x(2)=0.00001;
       end
LP=4.*(((29433./(G.*x)+(3300000.*b.^2)./(G.*x).^2).^0.5)./1000).^2;
[MSD,tau,weightval]=MSDcalc2(Spot(:,1)*pixelSize,Spot(:,2)*pixelSize,(Spot(:,9)-Spot(:,12))*frameTime);

%% Fit to MSD

[xData, yData, weights] = prepareCurveData( tau(1:nPoints-1), MSD(1:nPoints-1), weightval(1:nPoints-1) );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [0 LP(3)];
opts.Upper = [Inf LP(2)];
opts.Weights = weights;

% Fit model to data.
try
[fitresult, gof] = fit( xData, yData, ft, opts );
catch
    opts.Lower = [0 0];
opts.Upper = [Inf 100];
    [fitresult, gof] = fit( xData, yData, ft, opts );
end
fitCoeffs=coeffvalues(fitresult);
DiffusionConst=fitCoeffs(1)/4;
LocPrecision=(fitCoeffs(2).^0.5)./4;
% % Plot fit with data.
% figure( 'Name', 'untitled fit 2' );
% h = plot( fitresult, xData, yData );
% legend( h, 'MSD vs. tau with weightval', 'untitled fit 2', 'Location', 'NorthEast' );
% % Label axes
% xlabel tau
% ylabel MSD
% grid on

end