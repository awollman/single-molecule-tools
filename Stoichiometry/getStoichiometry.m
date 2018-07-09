function stoichiometry=getStoichiometry(spot, Isingle,noFrames,method,bleachTime)

                SpotT=spots(:,9)-mean(spots(:,12));
                SpotI=spots(:,5);
if noFrames==0
    noFramesToUse=size(spot,1);
end
switch method
    case 1 %Initial intensity
        stoichiometry=spot(1,5);
    case 2 %Mean intensity
        stoichiometry=mean(SpotI(1:noFramesToUse));
    case 3 %Linear intensity extrapolation
% determine stoichiometry by fitting a line and using intercept
pfit=polyfit(SpotT(1:noFramesToUse)-1,SpotI(1:noFramesToUse),1); %%REPLACE WITH CONSTRAINED FIT
    stoichiometry=pfit(2);
    case 4 %Exponential intensity extrapolation fixed bleach time
    case 5 %Exponential intensity extrapolation fitted bleach time
end
stoichiometry=stoichiometry/Isingle;
end