
function [trackArray,spotsInTracks]=trackAnalyser(spots,segmentedCell,tracksFile,cellNo,params)
%% function trackAnalyser
% does what you'd think, pulls any numbers from string in tracksFile and
% puts in column 7 of the array
%if PSFsize=0, no correction for PSF size or motion blur
pixelSize=params.pixelSize;
frameLimitS=params.frameLimitS;
frameLimitD=params.frameLimitD;
frameTime=params.frameTime;
PSFsize=params.PSFsize;
Isingle=params.Isingle
stoichMethod=params.stoichMethod;
bleachTime=params.bleachTime;
spotsInTracks=[];



cellCoord=[];

spotInd=ismember(round(spots(:,1:2)),cellCoord,'rows');
trajNo=unique(spots(spotInd,10));
trajNo(trajNo==0)=[];
for trajInd=1:length(trajNo)
    t=trajNo(trajInd);
    if min(spots(spots(:,10)==t,9))<(frameLimit+min(spots(:,9)))
        if length(spots(spots(:,10)==t,9))>max([frameLimitS, frameLimitD])
            trackNo=trackNo+1;
            SpotX=spots(spots(:,10)==t,1);
            SpotY=spots(spots(:,10)==t,2);
            SpotSX=spots(spots(:,10)==t,7);
            SpotSY=spots(spots(:,10)==t,6);
            trackArray(trackNo,1)=cellNo;
            
            [trackArray(trackNo,2),bleachTime]=getStoichiometry(spots(spots(:,10)==t,:), Isingle,frameLimitS,stoichMethod,bleachTime)
            %diffusion coefficient
            [trackArray(trackNo,3),MSD,tau,LocPrecision]=getDiffusion3(spots(spots(:,10)==t,:),frameTimeD,pixelSize,frameLimitD);
            % first displacement D(1)
            trackArray(trackNo,4)=((SpotX(2)-SpotX(1)).^2+(SpotY(2)-SpotY(1)).^2)*(pixelSize).^2/(4*frameTime);
            trackArray(trackNo,5)=(SpotSX(1)^2+SpotSY(1)^2)^0.5; %spot size
            trackArray(trackNo,6)=t; %trajectory number
            trackArray(trackNo,7)= str2num(tracksFile(regexp(tracksFile,'\d')));
            trackArray(trackNo,8:11)=0; % these will be assigned later with colocalisation info
            spotsInCell=spots(spots(:,10)==t,:);
            spotsInTracks=cat(1,spotsInTracks,spotsInCell);
        end
    end
end

end


