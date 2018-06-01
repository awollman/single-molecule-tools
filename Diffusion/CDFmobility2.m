function [fitresult1,fitresult2]=CDFmobility2(DiffusionConstNuc)
figure;
cdf=CDFeval(DiffusionConstNuc);
x=cdf(1,:);
y=cdf(2,:);
[xData, yData] = prepareCurveData( x, y );
[fitresult1, gof1,chi21] = expfit1(xData, yData);
[fitresult2, gof2,chi22] = expfit2(xData, yData);
[fitresult3, gof3,chi23] = expfit3(xData, yData);



subplot(1,3,1)
plot( fitresult1, xData, yData );
set(gca,'XScale','log');
%title(strcat('R2=',num2str(gof1.rsquare)))
title(strcat('reduced chi2=',num2str(chi21/(length(DiffusionConstNuc)-2))))
%xlabel(num2str(coeffvalues(fitresult1)))
[errorA1, errorD1]=errorCalc(fitresult1);
xlabel(figureLegend(fitresult1,errorA1, errorD1))

subplot(1,3,2)
plot( fitresult2, xData, yData );
set(gca,'XScale','log');
%title(strcat('R2=',num2str(gof2.rsquare)))
title(strcat('reduced chi2=',num2str(chi22/(length(DiffusionConstNuc)-4))))

%xlabel(num2str(coeffvalues(fitresult2)))
[errorA1, errorD1]=errorCalc(fitresult2);
xlabel(figureLegend(fitresult2,errorA1, errorD1))
subplot(1,3,3)
plot( fitresult1, xData, yData );
set(gca,'XScale','log');
%title(strcat('R2=',num2str(gof3.rsquare)))
title(strcat('reduced chi2=',num2str(chi23/(length(DiffusionConstNuc)-6))))
%xlabel(num2str(coeffvalues(fitresult3)))
[errorA1, errorD1]=errorCalc(fitresult3);
xlabel(figureLegend(fitresult3,errorA1, errorD1))
end

function [fitresult, gof,chi2] = expfit1(xData, yData)

% Set up fittype and options.
ft = fittype( 'A*(1-exp(-x/D))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.662381860399481 0.244165286790279];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
chi2=sum((yData(yData>0)-feval(fitresult,xData(yData>0))).^2./yData(yData>0).^2);
% Plot fit with data.
end

function [fitresult, gof,chi2] = expfit2(xData, yData)


% Set up fittype and options.
ft = fittype( 'A1*(1-exp(-x/D1))+A2*(1-exp(-x/D2))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 -Inf 0 0 0];
opts.Upper = [1 1 Inf Inf Inf Inf];
opts.StartPoint = [0.295507250831597 0.527846830418798 0.680178371230502 0.411593513407535];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
chi2=sum((yData(yData>0)-feval(fitresult,xData(yData>0))).^2./yData(yData>0).^2);
end

function [fitresult, gof,chi2] = expfit3(xData, yData)

ft = fittype( 'A1*(1-exp(-x/D1))+A2*(1-exp(-x/D2))+(A3)*(1-exp(-x/D3))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 -Inf 0 0 0];
opts.StartPoint = [0.295507250831597 0.527846830418798 0.602638218036397 0.680178371230502 0.411593513407535 0.750520055923736];
opts.Upper = [1 1 1 Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
chi2=sum((yData(yData>0)-feval(fitresult,xData(yData>0))).^2./yData(yData>0).^2);
end

function [errorA1, errorD1]=errorCalc(fitresult1)
coeffs1=coeffvalues(fitresult1);
confInts1=confint(fitresult1);
for n=1:length(coeffs1)/2
    errorA1(n)=abs((confInts1(1,n)-confInts1(2,n))/2);
    errorD1(n)=abs((confInts1(1,n)-confInts1(2,n))/2);
end
end

function [figureText]=figureLegend(fitresult1,errorA1, errorD1)
coeffs=coeffvalues(fitresult1);


n=1;
line1=strcat('A=',num2str(coeffs(n),2),setstr(177),...
    num2str(errorA1(n),2),' D=',num2str(coeffs(2*n),2),setstr(177),num2str(errorD1(n),2));
figureText=line1;
if length(coeffs)/2>1
n=2;
line2=strcat('A=',num2str(coeffs(n),2),setstr(177),...
    num2str(errorA1(n),2),' D=',num2str(coeffs(2*n),2),setstr(177),num2str(errorD1(n),2));
figureText={line1;line2};
end
if length(coeffs)/2>2
n=3;
line3=strcat('A=',num2str(coeffs(n),2),setstr(177),...
    num2str(errorA1(n),2),' D=',num2str(coeffs(2*n),2),setstr(177),num2str(errorD1(n),2));
figureText={line1;line2;line3};
end

end