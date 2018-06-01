%% Initialise
noSpots=1;
Isingle=4000;
BGmean=500; %mean background pixel intensity
BGstd=120; %standard deviation of background pixels
noFrames=5000;
sizeN=50;
sizeM=50;
bleachTime=10; %in frames, if 0 then no bleaching
diffusionCoeff=0;%um2/s
nDiffPoints=4; %number of MSD points to calculate diffusion const
frameTime=0.005; %seconds
pixelSize=0.120; %um
PSFwidth=0.160/pixelSize; %sigma of a Gaussian, ~2/3 airy disk diameter
fileName='test.tif';
%make a spot array the same size as normal
spotsReal=zeros(noSpots*noFrames,12);
spotsReal(:,5)=Isingle;

% initialise the spot co-ords
spotsReal(1:noSpots,1)=rand(noSpots,1)*sizeN;
spotsReal(1:noSpots,2)=rand(noSpots,1)*sizeM;
spotsReal(1:noSpots,9)=1;
spotsReal(1:noSpots,10)=1:noSpots;

%% simulate diffusion
currentTrajNo=1:noSpots;
S = sqrt(2*diffusionCoeff*frameTime)/pixelSize;
for t=1:noFrames-1
    spotsReal(noSpots*t+1:noSpots*(t+1),9)=t+1;
    spotsReal(noSpots*t+1:noSpots*(t+1),10)=currentTrajNo;
    spotsReal(noSpots*t+1:noSpots*(t+1),1)=normrnd(spotsReal(noSpots*(t-1)+1:noSpots*t,1),S*ones(noSpots,1),noSpots,1);
    spotsReal(noSpots*t+1:noSpots*(t+1),2)=normrnd(spotsReal(noSpots*(t-1)+1:noSpots*t,2),S*ones(noSpots,1),noSpots,1);
    
    if bleachTime>0
        bleachedSpots=currentTrajNo(rand(1,noSpots)<1/bleachTime);      
        for b=1:length(bleachedSpots)
            spotsReal(spotsReal(:,9)==t+1 & spotsReal(:,10)==bleachedSpots(b),1)=rand(1)*sizeN;
            spotsReal(spotsReal(:,9)==t+1 & spotsReal(:,10)==bleachedSpots(b),2)=rand(1)*sizeM;
            spotsReal(spotsReal(:,9)==t+1 & spotsReal(:,10)==bleachedSpots(b),10)=max(spotsReal(:,10))+1;            
        end
    end
    
        currentTrajNo=spotsReal(noSpots*t+1:noSpots*(t+1),10);
end

%% simulate an image stack and save
[Xpos,Ypos] = meshgrid(1:sizeN,1:sizeM);
GaussFrame=zeros(sizeN,sizeM,noFrames);
for t=1:noFrames
    for spts=1:noSpots
        spotInd=(t-1)*noSpots+spts;
        GaussFrame(:,:,t)=GaussFrame(:,:,t)+(spotsReal(spotInd,5)./(2.*pi.*PSFwidth.^2))*...
            exp(-(((Xpos-spotsReal(spotInd,1)).^2)./(2.*PSFwidth^2)+((Ypos-spotsReal(spotInd,2)).^2)./(2.*PSFwidth^2)));
    end
    %add poisson noise to simulate shot noise on the fluorescent spots
    GaussFrame(:,:,t)=imnoise(uint16(GaussFrame(:,:,t)),'poisson');
    %add Gaussian noise to simulate background camera noise
    GaussFrame(:,:,t)=GaussFrame(:,:,t)+normrnd(BGmean,BGstd,size(GaussFrame(:,:,t)));
    if t==1
        imwrite(uint16(GaussFrame(:,:,t)),fileName)
    else
        imwrite(uint16(GaussFrame(:,:,t)),fileName,'WriteMode','append')
    end
end
%% Track
createP
p.show_output=0;

[SpotsCh1, ~, ~,~, ~, ~,~] = ADEMScode2_84(GaussFrame,p);

for t=1:max(SpotsCh1(:,10))
    if sum(SpotsCh1(:,10)==t)>=nDiffPoints
        [DiffusionConst(t),MSD,tau,LocPrecisionDiff(t)]=getDiffusion3(SpotsCh1(SpotsCh1(:,10)==t,:),frameTime,pixelSize,nDiffPoints);
    end
end

for f=1:noFrames
    trackSpotsFrame=SpotsCh1(SpotsCh1(:,9)==f,:);
    realSpotsFrame=spotsReal(spotsReal(:,9)==f,:);
    numErroneousSpots(f)=0;
    for k=1:size(trackSpotsFrame,1)
        if min(((realSpotsFrame(:,1)-trackSpotsFrame(k,1)).^2+(realSpotsFrame(:,2)-trackSpotsFrame(k,2)).^2).^0.5)>2
            trackSpotsFrame(k,10)=-1;
            numErroneousSpots(f)=numErroneousSpots(f)+1;
        end
    end
    numSpotsMissed(f)=0;
    for j=1:noSpots
        if min((trackSpotsFrame(:,1)-realSpotsFrame(j,1)).^2+(trackSpotsFrame(:,2)-realSpotsFrame(j,2)).^2)>2
            
            numSpotsMissed(f)=numSpotsMissed(f)+1;
        end
        
    end
end


%% plot things
figure;
subplot(2,2,1)
imshow(GaussFrame(:,:,end),[])
hold on
for s=1:max(spotsReal(:,10)),plot(spotsReal(spotsReal(:,10)==s,1),spotsReal(spotsReal(:,10)==s,2)), end
title('last frame with real tracks')

subplot(2,2,2)
imshow(GaussFrame(:,:,end),[])
hold on
for s=1:max(SpotsCh1(:,10)),plot(SpotsCh1(SpotsCh1(:,10)==s,1),SpotsCh1(SpotsCh1(:,10)==s,2)), end
title('last frame with tracked tracks')

subplot(2,2,3)
[counts,x]=KDFplot(DiffusionConst);
hold on
plot(ones(length(counts)/10)*diffusionCoeff,counts(1:10:end),'--k')
legend('Measured diffusion coeffs','Simulated diffusion coeff')
xlabel('Diffusion Coeff {\mum}^2/s')
ylabel('probability')

subplot(2,2,4)
bar([mean(numSpotsMissed),mean(numErroneousSpots)])
hold on
errorbar(1:2,[mean(numSpotsMissed),mean(numErroneousSpots)],[std(numSpotsMissed),std(numErroneousSpots)],'k')
xticklabels({'false negartives','false positives'})
ylabel('Number of spots')