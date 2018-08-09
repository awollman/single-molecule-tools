% links spots and calculates overlap integral
% links spots in Ch2 to spots in Ch1
% adds columns to spot array:
%column 13=index of linked spot in other spot array
%column 14=overlap integral for that spot
%column 15=trajectory number of linked spot in other channel

function [SpotsCh1linked, SpotsCh2linked]=Colocaliser2(SpotsCh1,SpotsCh2,params)

%% define parameters
% params.transform=[0,0];
if nargin<3
    % any co-ordinate transformation between the two channels
    params.transform=[0,0];
    % maximum distance in pixels to link spots
    params.d=5;
    % maximum overlap integral to link spots
    params.overlap=0.75;
    
    % set if=0 link all spots regardless of frame
    %     if=1 link only spots in the same frame
    %     if=2 link spots in alternate frames (for ALEX)
    params.frameLinkMethod=0;
    
    % set so that spots are only assigned 1 partner
    params.Unique=1;
    
    %set to display graphs
    params.showOutput=0;
end

%% perform co-ordinate transforms between channels here

SpotsCh2(:,1)=SpotsCh2(:,1)+params.transform(1);
SpotsCh2(:,2)=SpotsCh2(:,2)+params.transform(2);

%% link spots
switch params.frameLinkMethod
    case 0
        Spots1=SpotsCh1;
        Spots2=SpotsCh2;
        [SpotsCh1linked, SpotsCh2linked]=linker(Spots1, Spots2);
    case 1 %link spots in same frame
        SpotsCh1linked=[];
        SpotsCh2linked=[];
        for x=min([SpotsCh1(:,9)',SpotsCh2(:,9)']):max([SpotsCh1(:,9)',SpotsCh2(:,9)'])
            Spots1=SpotsCh1(SpotsCh1(:,9)==x,:);
            Spots2=SpotsCh2(SpotsCh2(:,9)==x,:);
            [SpotsCh1linkedTEMP, SpotsCh2linkedTEMP]=linker(Spots1, Spots2);
            % NEED TO ADD TO COL13 SO SPOT NO CONSISTENT
            if isempty(SpotsCh1linked)==0
                SpotsCh1linkedTEMP(SpotsCh1linkedTEMP(:,13)>0,13)=SpotsCh1linkedTEMP(SpotsCh1linkedTEMP(:,13)>0,13)+size(SpotsCh2linked,1);
                SpotsCh2linkedTEMP(SpotsCh2linkedTEMP(:,13)>0,13)=SpotsCh2linkedTEMP(SpotsCh2linkedTEMP(:,13)>0,13)+size(SpotsCh1linked,1);
            end
            
            SpotsCh1linked=cat(1,SpotsCh1linked,SpotsCh1linkedTEMP);
            SpotsCh2linked=cat(1,SpotsCh2linked,SpotsCh2linkedTEMP);
        end
    case 2 %link alternate frames ALEX
          SpotsCh1linked=[];
        SpotsCh2linked=[];
        [firstFrame, firstFrameChannel]=min([min(SpotsCh1(:,9)),min(SpotsCh2(:,9))]);
        [lastFrame, lastFrameChannel]=max([max(SpotsCh1(:,9)),max(SpotsCh2(:,9))]);
        for x=firstFrame:2:lastFrame
            if firstFrameChannel==1
                frame1=x;
                frame2=x+1;
            else
                frame1=x+1;
                frame2=x;
            end
            Spots1=SpotsCh1(SpotsCh1(:,9)==frame1,:);
            Spots2=SpotsCh2(SpotsCh2(:,9)==frame2,:);
            [SpotsCh1linkedTEMP, SpotsCh2linkedTEMP]=linker(Spots1, Spots2);
            % NEED TO ADD TO COL13 SO SPOT NO CONSISTENT
            if isempty(SpotsCh1linked)==0
                SpotsCh1linkedTEMP(SpotsCh1linkedTEMP(:,13)>0,13)=SpotsCh1linkedTEMP(SpotsCh1linkedTEMP(:,13)>0,13)+size(SpotsCh2linked,1);
                SpotsCh2linkedTEMP(SpotsCh2linkedTEMP(:,13)>0,13)=SpotsCh2linkedTEMP(SpotsCh2linkedTEMP(:,13)>0,13)+size(SpotsCh1linked,1);
            end
            
            SpotsCh1linked=cat(1,SpotsCh1linked,SpotsCh1linkedTEMP);
            SpotsCh2linked=cat(1,SpotsCh2linked,SpotsCh2linkedTEMP);
        end
end


if params.showOutput==1
    try
    %% scatter plot
    figure
    subplot(1,3,1)
    scatter(SpotsCh1linked(SpotsCh1linked(:,14)>params.overlap,1),SpotsCh1linked(SpotsCh1linked(:,14)>params.overlap,2),'g')
    hold on
    scatter(SpotsCh2linked(SpotsCh2linked(:,14)>params.overlap,1),SpotsCh2linked(SpotsCh2linked(:,14)>params.overlap,2),'r')
    xlabel('n (pixels)')
    ylabel('m (pixels)')
    title('scatter co-localised co-ordinates')
    xlim([0,round(max(SpotsCh1(:,1)))])
    ylim([0,round(max(SpotsCh1(:,2)))])
   
    
    %% Distribution of overlap integrals
    subplot(1,3,2)
    KDFplot(SpotsCh1linked(SpotsCh1linked(:,14)>params.overlap,14));
    xlabel('Overlap Integrals > params.overlap')
    
    %% reconstructed image
    ScaleFactor=1;
    GaussFrame=zeros(ScaleFactor*[round(max(SpotsCh1(:,2))),round(max(SpotsCh1(:,1)))]);
    GaussFrame2=GaussFrame;
    [Xpos,Ypos] = meshgrid(1:size(GaussFrame,2),1:size(GaussFrame,1));
    
    spotInd=find(SpotsCh1linked(:,14)>params.overlap);
    for p=1:length(find(SpotsCh1linked(:,14)>params.overlap))
        t=spotInd(p);
        Intensity_rec=SpotsCh1linked(t,5);
        GaussFrame=GaussFrame+(Intensity_rec./(2.*pi.*SpotsCh1linked(t,6).*SpotsCh1linked(t,7)))*exp(-(((Xpos-ScaleFactor*SpotsCh1linked(t,1)).^2)./(2.*SpotsCh1linked(t,6)^2)+((Ypos-ScaleFactor*SpotsCh1linked(t,2)).^2)./(2.*SpotsCh1linked(t,7)^2)));
    end
    
    spotInd=find(SpotsCh2linked(:,14)>params.overlap);
    for p=1:length(find(SpotsCh2linked(:,14)>params.overlap))
        t=spotInd(p);
        Intensity_rec=SpotsCh2linked(t,5);
        GaussFrame2=GaussFrame2+(Intensity_rec./(2.*pi.*SpotsCh2linked(t,6).*SpotsCh2linked(t,7)))*exp(-(((Xpos-ScaleFactor*SpotsCh2linked(t,1)).^2)./(2.*SpotsCh2linked(t,6)^2)+((Ypos-ScaleFactor*SpotsCh2linked(t,2)).^2)./(2.*SpotsCh2linked(t,7)^2)));
    end
    subplot(1,3,3)
    imshow(cat(3,mat2gray(GaussFrame2),mat2gray(GaussFrame), zeros(size(GaussFrame))),[])
    title('re-constructed co-localised image')
    catch
        disp('colocalisation plotting failed, likely because there are no co-linked spots')
    end
   % pause
end

% does the business of linking spots
    function [SpotsCh1linked, SpotsCh2linked]=linker(Spots1, Spots2)
        SpotsCh1linked=Spots1;
        SpotsCh2linked=Spots2;
        distanceMatrix=pdist2([Spots1(:,1),Spots1(:,2)],[Spots2(:,1),Spots2(:,2)]);
        [Spot1ind,Spot2ind] = find(distanceMatrix<params.d); %put params.d here
        if isempty(Spot1ind)==0
        sigma1sq=(Spots1(Spot1ind,6).^2+Spots1(Spot1ind,7).^2);
        sigma2sq=(Spots2(Spot2ind,6).^2+Spots2(Spot2ind,7).^2);
         overlapIntMatrix=zeros(size(distanceMatrix));
        % vector containing all overlap integrals
        try %this prevents error at low spot numbers where vectors are longer in rows/cols
            %seems to be a weird matlab feature which switches them around
            overlapInt=exp(-distanceMatrix(distanceMatrix<params.d).^2./(2.*(sigma1sq(:)+sigma2sq(:))));
            overlapIntMatrix(distanceMatrix<params.d)=overlapInt;
        catch
            overlapInt=exp(-distanceMatrix(distanceMatrix<params.d)'.^2./(2.*(sigma1sq+sigma2sq)));
            overlapIntMatrix(distanceMatrix<params.d)=overlapInt;
        end
       
        % matix same size as distance matrix with overlap integral at each
        % distance <5
   %     overlapIntMatrix(distanceMatrix<params.d)=overlapInt;
        [maxValCh1, maxICh1]=max(overlapIntMatrix,[],1);
        
        maxICh2=1:length(maxICh1);
        maxICh1(maxValCh1< params.overlap)=[];
        maxICh2(maxValCh1< params.overlap)=[];
        maxValCh1(maxValCh1< params.overlap)=[];
        % choose the highest overlap integral for non-unique assigments
        if  params.Unique==1
            C=unique(maxICh1);
            for s=1:length(unique(maxICh1))
                if sum(maxICh1==C(s))>1
                    deleteIndex=maxICh1==C(s) & maxValCh1~=max(maxValCh1(maxICh1==C(s)));
                    maxICh1(deleteIndex)=[];
                    maxValCh1(deleteIndex)=[];
                    maxICh2(deleteIndex)=[];
                end
            end
        end
        % assign the overlap integral and spot indices to spot arrays
        SpotsCh2linked(maxICh2,13)=maxICh1;
        SpotsCh2linked(maxICh2,14)=maxValCh1;
        SpotsCh1linked(maxICh1,13)=maxICh2;
             

        SpotsCh1linked(maxICh1,14)=maxValCh1;
        SpotsCh1linked(maxICh1,15)=SpotsCh2linked(maxICh2,10);
        SpotsCh2linked(maxICh2,15)=SpotsCh1linked(maxICh1,10);
        
        SpotsCh1linked(maxICh1,16)=((SpotsCh1linked(maxICh1,1)-SpotsCh2linked(maxICh2,1)).^2+(SpotsCh1linked(maxICh1,2)-SpotsCh2linked(maxICh2,2)).^2).^0.5;
        SpotsCh2linked(maxICh2,16)=((SpotsCh1linked(maxICh1,1)-SpotsCh2linked(maxICh2,1)).^2+(SpotsCh1linked(maxICh1,2)-SpotsCh2linked(maxICh2,2)).^2).^0.5;
        else
        SpotsCh1linked(:,13:16)=0;
        SpotsCh2linked(:,13:16)=0;

        end
    end



end