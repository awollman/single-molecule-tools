
function [trackArray,spotsInTracks]=trackAnalyser(spots,segmentation,tracksFile,cellNo,params)
%% function trackAnalyser
% does what you'd think, pulls any numbers from string in tracksFile and
% puts in column 7 of the array


if nargin<5
    params.pixelSize=0.08; %pixel size in microns
    params.frameLimitS=4; %num frames to use in stoich
    params.frameLimitD=4; %num frames to use in diffusion
    params.frameTime=0.005; % time between frames in seconds
    params.Isingle=5000; %characteristic intensity of a single fluorophore
    params.stoichMethod=1; %method for calculating stoichiomtry, see getStoichiometry
    params.bleachTime=5; %required for some stoich methods
    params.showOutput=1; %makes plots of each cell
    params.frameLimitAll=10; %number of frames in include in analysis
end

pixelSize=params.pixelSize;
frameLimitS=params.frameLimitS;
frameLimitD=params.frameLimitD;
frameTime=params.frameTime;
Isingle=params.Isingle;
stoichMethod=params.stoichMethod;
bleachTime=params.bleachTime;
spotsInTracks=[];



spots(spots(:,10)==0,:)=[];
trackArray=[];
if params.showOutput==1
    figure;
    subplot(2,3,1)
    imshow(sum(segmentation,3),[])
    hold on
end

if size(segmentation,3)==1
    
segmentedCell=segmentation;
else
    segmentedCell=segmentation(:,:,1);
for s=1:size(segmentation,3)-1
            rp=regionprops(segmentation(:,:,s+1),'centroid');
        text(rp.Centroid(1), rp.Centroid(2),num2str(s),'color','red')
end

[cellCoord(:,2), cellCoord(:,1)]=find(segmentedCell);
spotInd=ismember(round(spots(:,1:2)),cellCoord,'rows');
trajNo=unique(spots(spotInd,10));
trajNo(trajNo==0)=[];
trackNo=0;
for trajInd=1:length(trajNo)
    t=trajNo(trajInd);
    if min(spots(spots(:,10)==t,9))<(params.frameLimitAll+min(spots(:,9)))
        if length(spots(spots(:,10)==t,9))>max([frameLimitS, frameLimitD])
            trackNo=trackNo+1;
            SpotX=spots(spots(:,10)==t,1);
            SpotY=spots(spots(:,10)==t,2);
            SpotSX=spots(spots(:,10)==t,7);
            SpotSY=spots(spots(:,10)==t,6);
            trackArray(trackNo,1)=cellNo;
            %stoichiometry
            [trackArray(trackNo,2),bleachTime]=getStoichiometry(spots(spots(:,10)==t,:), Isingle,frameLimitS,stoichMethod,bleachTime);
            %diffusion coefficient
            [trackArray(trackNo,3),MSD,tau,LocPrecision]=getDiffusion3(spots(spots(:,10)==t,:),frameTime,pixelSize,frameLimitD);
            % first displacement D(1)
            trackArray(trackNo,4)=((SpotX(2)-SpotX(1)).^2+(SpotY(2)-SpotY(1)).^2)*(pixelSize).^2/(4*frameTime);
            trackArray(trackNo,5)=(SpotSX(1)^2+SpotSY(1)^2)^0.5; %spot size
            trackArray(trackNo,6)=t; %trajectory number
            trackArray(trackNo,7)= str2num(tracksFile(regexp(tracksFile,'\d')));
            trackArray(trackNo,8:11)=0; % these will be assigned later with colocalisation info
            trackArray(trackNo,12)=min(spots(spots(:,10)==t,9)); %first frame of track
            trackArray(trackNo,13)=max(spots(spots(:,10)==t,9)); %last frame of track
                trackArray(trackNo,14)=0;
                        trackArray(trackNo,15)=0;
            if size(segmentation,3)==1

            else
               
                for s=1:size(segmentation,3)-1
                    clear compCoord
                    [compCoord(:,2), compCoord(:,1)]=find(segmentation(:,:,s+1));
                    spotInd2=ismember(round(spots(spots(:,10)==t,1:2)),compCoord,'rows');
                    if sum(spotInd2)>0
                        trackArray(trackNo,14)=s;
                        trackArray(trackNo,15)=sum(spotInd2)/sum(spots(:,10)==t);
                        if sum(spotInd2(1:floor(end/2)))/sum(spotInd2(ceil(end/2):end))
                            % positive if starts in compartment
                        trackArray(trackNo,15)=trackArray(trackNo,15)*1;
                        else
                            % negative if starts outside compartment
                        trackArray(trackNo,15)=trackArray(trackNo,15)*-1;
                        end
                    end
                end
            end
            spotsInCell=spots(spots(:,10)==t,:);
            spotsInTracks=cat(1,spotsInTracks,spotsInCell);
            
            if params.showOutput==1
                
                subplot(2,3,1)
                plot(SpotX,SpotY)
                title(strcat('Cell number= ',num2str(cellNo)))
                subplot(2,3,2)
                plot(spots(spots(:,10)==t,9),spots(spots(:,10)==t,5)/params.Isingle)
                hold on
                xlabel('frame number')
                ylabel('Intensity (number of fluorophores)')
            end
        end
    end
end

if params.showOutput==1
    try
    subplot(2,3,3)
[counts,x]=ksdensity(trackArray(:,2));
plot(x,counts)
xlabel('Stoichiometry')
ylabel('Probability')
    subplot(2,3,4)
[counts,x]=ksdensity(trackArray(:,3));
plot(x,counts)
xlabel('Diffusion Coeff')
ylabel('Probability')
subplot(2,3,5)
scatter(trackArray(:,2),trackArray(:,3),'.')
xlabel('Stoichiometry')
ylabel('Diffusion Coeff')
subplot(2,3,6)
[counts,x]=histcounts(spots(:,10),'binwidth',1);
histogram(counts(counts>0),'binwidth',1)
hold on
[counts,x]=histcounts(spotsInTracks(:,10),'binwidth',1);
histogram(counts(counts>0),'binwidth',1)
xlabel('Track length (frames)')
ylabel('Frequency')
legend('All spots','Spots in cell and frame limit')
    catch
        disp('Plotting error')
    end
end


end


