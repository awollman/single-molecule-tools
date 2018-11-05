function [stoichiometry,bleachTime]=getStoichiometry(spot, Isingle,noFramesToUse,method,bleachTime)

SpotT=spot(:,9)-mean(spot(:,12));
SpotI=spot(:,5);
if noFramesToUse==0
    noFramesToUse=size(spot,1);
end
switch method
    case 1 %Initial intensity
        stoichiometry=spot(1,5);
    case 2 %Mean intensity
        stoichiometry=mean(SpotI(1:noFramesToUse));
    case 3 %Linear intensity extrapolation
        % determine stoichiometry by fitting a line and using intercept
        [xData, yData] = prepareCurveData( SpotT(1:noFramesToUse), SpotI(1:noFramesToUse) );
        % Set up fittype and options.
        ft = fittype( 'poly1' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Lower = [-Inf spot(1,5)/2];
        opts.Upper = [0 Inf];
        
        [fitresult, gof] = fit( xData, yData, ft, opts );
        
        fitCoeffs=coeffvalues(fitresult);
        stoichiometry=fitCoeffs(2);
        bleachTime=fitCoeffs(1);
        %         pfit=polyfit(SpotT(1:noFramesToUse)-1,SpotI(1:noFramesToUse),1); %%REPLACE WITH CONSTRAINED FIT
        %         stoichiometry=pfit(2);
    case 4 %Exponential intensity extrapolation fixed bleach time
        stoichiometry=SpotI(1:noFramesToUse)\exp(-SpotT(1:noFramesToUse)/bleachTime);
    case 5 %Exponential intensity extrapolation fitted bleach time
        [xData, yData] = prepareCurveData( SpotT(1:noFramesToUse), SpotI(1:noFramesToUse) );
        
        % Set up fittype and options.
        ft = fittype( 'I0*exp(-b*x)+Is', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [5*Isingle Isingle 0.1];
        opts.Lower = [0 0 0];
        opts.Upper = [Inf Isingle*2 Inf];
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
        coeffvals = coeffvalues(fitresult);
        stoichiometry=coeffvals(1);
        bleachTime=coeffvals(3);
end
stoichiometry=stoichiometry/Isingle;
end