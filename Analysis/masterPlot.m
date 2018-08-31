function h1=masterPlot(trackArray,params,plotSelect,h1)

    %select to do various extra analysis
    %=0 plots nothing!
    %=1 normal Stoichiometry and diffusion etc..
    %=2 colocalisation analysis different plots
    %=3 colocalisation analysis same plots
    %=5 periodicity analysis
    %=6 CDF analysis
    
    if ~exist('params.IgnoreEmpty')
        params.IgnoreEmpty=0;
    end


cmap=colormap('jet');
maxColor=5;
            plotName='Track Plot';
            if exist('h1')
                figure(h1)
                hchild1=get(h1,'children');
                axisProps=findobj(hchild1(end),'Type','line');
                colorInd=size(axisProps,1)+1;
            else
                colorInd=1;
                if plotSelect==0
                    h1=[];
                else
                    h1=figure('Name',plotName,'NumberTitle','off');
                end
            end
            Isingle=params.Isingle;

    %% means
    disp(strcat('Mean stoichiometry=',num2str(mean(trackArray(:,2)),2),'+/-',num2str(std(trackArray(:,2)),2)))
    disp(strcat('Mean D=',num2str(mean(trackArray(:,3)),2),'+/-',num2str(std(trackArray(:,3)),2)))
%     disp(strcat('N unlinked trajectories=',num2str(sum(trackArray(:,8)==0))))
%     disp(strcat('N linked trajectories=',num2str(sum(trackArray(:,8)>0))))
    if  params.IgnoreEmpty==0
    for c=1:max(trackArray(:,1))
         unlinkTrack(c)=sum(trackArray(trackArray(:,1)==c,8)==0);
        linkTrack(c)=sum(trackArray(trackArray(:,1)==c,8)>0);

    end
    disp(strcat('N cells=',num2str(max(trackArray(:,1)))))
    else
        uniqueCellNo=unique([trackArrayCh1(:,1)',trackArrayCh2(:,1)']);
     for c=1:length(uniqueCellNo)
         unlinkTrack(c)=sum(trackArray(trackArray(:,1)==uniqueCellNo(c),8)==0);
        linkTrack(c)=sum(trackArray(trackArray(:,1)==uniqueCellNo(c),8)>0);

     end
    disp(strcat('N cells=',num2str(length(uniqueCellNo))))
    end
     disp(strcat('N unlinked trajectories/cell=',num2str(mean(unlinkTrack)),'+/-',num2str(std(unlinkTrack)/c^0.5)))
    disp(strcat('N linked trajectories/cell=',num2str(mean(linkTrack)),'+/-',num2str(std(linkTrack)/c^0.5)))

   disp(strcat('Mean unlinked stoichiometry=',num2str(mean(trackArray(trackArray(:,8)==0,2)),2),'+/-',num2str(std(trackArray(trackArray(:,8)==0,2)),2)))
            disp(strcat('Mean linked stoichiometry=',num2str(mean(trackArray(trackArray(:,8)>0,2)),2),'+/-',num2str(std(trackArray(trackArray(:,8)>0,2)),2)))
            disp(strcat('Mean unlinked D=',num2str(mean(trackArray(trackArray(:,8)==0,3)),2),'+/-',num2str(std(trackArray(trackArray(:,8)==0,3)),2)))
            disp(strcat('Mean linked D=',num2str(mean(trackArray(trackArray(:,8)>0,3)),2),'+/-',num2str(std(trackArray(trackArray(:,8)>0,3)),2)))
            disp(strcat('Mean linked seperation=',num2str(mean(trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000),3),'+/-',num2str(std(trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000),2)))
  
            %     subplot(2,4,8)
    %     str=strcat('Channel=',num2str(ch));
    %     text(0.5,0.5,str); axis off
    %annotation('textbox','String',str);
    
    
    %% more output
    
    %% Hardcore Analysis
    switch plotSelect
        case 0

        case 1
            
            
            
            %% plots
            subplot(2,4,1)
            %figure('Name','Stoichiometry Distribution','NumberTitle','off')
            [plotH,~,x]=KDFplotHlocal(trackArray(:,2),0.7);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            hold on
            xlim([min(x),max(x)])
            
            subplot(2,4,2)
            %figure('Name','Stoichiometry Histogram','NumberTitle','off')
            histogram(trackArray(:,2),0.5:1:round(max(trackArray(:,2)))+0.5,'FaceColor',cmap(round(64*(colorInd)/maxColor),:));
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            hold on
            %   xlim([0,150])
            subplot(2,4,3)
            %figure('Name','Diffusion Const Distribution','NumberTitle','off')
            [plotH,~,~]=KDFplotHlocal(trackArray(:,3),0.05);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            
            xlabel('Diffusion Const (um2/s)')
            ylabel('Probability Density')
            %xlim([-0.2,0.6])
            hold on
            subplot(2,4,4)
            %figure('Name','S vs D','NumberTitle','off')
            %scatter(AllDiff, AllmeanI);
            scatter(trackArray(:,3),trackArray(:,2),50,cmap(round(64*(colorInd)/maxColor),:),'filled');
            %xlabel('Diffusion Const (um2/s)')
            xlabel('Diffusion Const (um2/s)')
            ylabel('Stoichiometry')
            %xlim([0,0.6])
            hold on
            
            %no spots/cell
            subplot(2,4,5)
            [counts,x]=histcounts(trackArray(:,1));
            histogram(counts,'FaceColor',cmap(round(64*(colorInd)/maxColor),:));
            xlabel('Number of trajectories/cell')
            ylabel('Frequency')
            hold on
            %total number of molecules in spots/cell
            subplot(2,4,6)
            for c=1:max(trackArray(:,1))
                molSpotPerCell(c)=sum(trackArray(trackArray(:,1)==c,2));
            end
            [plotH,~,~]=KDFplotHlocal(molSpotPerCell);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            
            xlabel('total number of molecules in trajectories/cell')
            ylabel('Probability Density')
            hold on
            %mean spot stoich/cell
            subplot(2,4,7)
            for c=1:max(trackArray(:,1))
                molSpotPerCell(c)=mean(trackArray(trackArray(:,1)==c,2));
            end
            [plotH,~,~]=KDFplotHlocal(molSpotPerCell);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            
            xlabel('mean number of molecules in trajectory/cell')
            ylabel('Probability Density')
            hold on
        case 2
            %% plots
            disp(strcat('Mean unlinked stoichiometry=',num2str(mean(trackArray(trackArray(:,8)==0,2)),2),'+/-',num2str(std(trackArray(trackArray(:,8)==0,2)),2)))
            disp(strcat('Mean linked stoichiometry=',num2str(mean(trackArray(trackArray(:,8)>0,2)),2),'+/-',num2str(std(trackArray(trackArray(:,8)>0,2)),2)))
            disp(strcat('Mean unlinked D=',num2str(mean(trackArray(trackArray(:,8)==0,3)),2),'+/-',num2str(std(trackArray(trackArray(:,8)==0,3)),2)))
            disp(strcat('Mean linked D=',num2str(mean(trackArray(trackArray(:,8)>0,3)),2),'+/-',num2str(std(trackArray(trackArray(:,8)>0,3)),2)))
            disp(strcat('Mean linked seperation=',num2str(mean(trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000),3),'+/-',num2str(std(trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000),2)))
            subplot(2,4,1)
            %figure('Name','Stoichiometry Distribution','NumberTitle','off')
            [plotH,~,x]=KDFplotHlocal(trackArray(trackArray(:,8)==0,2),0.7);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            title('Unlinked Trajectories')
            hold on
            xlim([min(x),max(x)])
            subplot(2,4,2)
            %figure('Name','Stoichiometry Distribution','NumberTitle','off')
            [plotH,~,x]=KDFplotHlocal(trackArray(trackArray(:,8)>0,2),0.7);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            title('Linked Trajectories')
            hold on
            xlim([min(x),max(x)])
            subplot(2,4,3)
            %figure('Name','Diffusion Const Distribution','NumberTitle','off')
            [plotH,~,~]=KDFplotHlocal(trackArray(trackArray(:,8)==0,3),0.05);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            xlabel('Diffusion Const (um2/s)')
            ylabel('Probability Density')
            title('Unlinked Trajectories')
            %xlim([-0.2,0.6])
            hold on
            subplot(2,4,4)
            %figure('Name','Diffusion Const Distribution','NumberTitle','off')
            [plotH,~,~]=KDFplotHlocal(trackArray(trackArray(:,8)>0,3),0.05);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            xlabel('Diffusion Const (um2/s)')
            ylabel('Probability Density')
            title('linked Trajectories')
            %xlim([-0.2,0.6])
            hold on
            subplot(2,4,5)
            [counts,x]=histcounts(trackArray(trackArray(:,8)==0,1));
            histogram(counts,'FaceColor',cmap(round(64*(colorInd)/maxColor),:));
            xlabel('Number of unlinked trajectories/cell')
            ylabel('Frequency')
            hold on
            subplot(2,4,6)
            [counts,x]=histcounts(trackArray(trackArray(:,8)>0,1));
            histogram(counts,'FaceColor',cmap(round(64*(colorInd)/maxColor),:));
            xlabel('Number of linked trajectories/cell')
            ylabel('Frequency')
            hold on
            subplot(2,4,7)
            [plotH,~,~]=KDFplotHlocal(trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000,10);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            xlabel('Seperation between linked spots (nm)')
            ylabel('Probability Density')
            title('linked Trajectories')
            hold on
            subplot(2,4,8)
            title('linked Trajectories')
            scatter(trackArray(trackArray(:,11)>0,2),trackArray(trackArray(:,11)>0,11)*params.pixelSize*1000,50,cmap(round(64*(colorInd)/maxColor),:),'filled');
            hold on
            ylabel('Seperation between linked spots (nm)')
            xlabel('Stoichiometry')
        case 3
            %% plots
            disp(strcat('Mean unlinked stoichiometry=',num2str(mean(trackArray(trackArray(:,8)==0,2)),2),'+/-',num2str(std(trackArray(trackArray(:,8)==0,2)),2)))
            disp(strcat('Mean linked stoichiometry=',num2str(mean(trackArray(trackArray(:,8)>0,2)),2),'+/-',num2str(std(trackArray(trackArray(:,8)>0,2)),2)))
            disp(strcat('Mean unlinked D=',num2str(mean(trackArray(trackArray(:,8)==0,3)),2),'+/-',num2str(std(trackArray(trackArray(:,8)==0,3)),2)))
            disp(strcat('Mean linked D=',num2str(mean(trackArray(trackArray(:,8)>0,3)),2),'+/-',num2str(std(trackArray(trackArray(:,8)>0,3)),2)))
            
            subplot(2,2,1)
            %figure('Name','Stoichiometry Distribution','NumberTitle','off')
            [plotH,~,x]=KDFplotHlocal(trackArray(trackArray(:,8)==0,2),0.7);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            title('Unlinked Trajectories')
            hold on
            xlim([min(x),max(x)])
            % subplot(2,4,2)
            %figure('Name','Stoichiometry Distribution','NumberTitle','off')
            [plotH,~,x]=KDFplotHlocal(trackArray(trackArray(:,8)>0,2),0.7);
            set(plotH,'Color',cmap(round(64*(colorInd+1)/maxColor),:))
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            title('Linked Trajectories')
            hold on
            xlim([min(x),max(x)])
            subplot(2,2,2)
            %figure('Name','Diffusion Const Distribution','NumberTitle','off')
            [plotH,~,~]=KDFplotHlocal(trackArray(trackArray(:,8)==0,3),0.05);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            xlabel('Diffusion Const (um2/s)')
            ylabel('Probability Density')
            title('Unlinked Trajectories')
            %xlim([-0.2,0.6])
            hold on
            %  subplot(2,4,4)
            %figure('Name','Diffusion Const Distribution','NumberTitle','off')
            [plotH,~,~]=KDFplotHlocal(trackArray(trackArray(:,8)>0,3),0.05);
            set(plotH,'Color',cmap(round(64*(colorInd+1)/maxColor),:))
            xlabel('Diffusion Const (um2/s)')
            ylabel('Probability Density')
            title('linked Trajectories')
            %xlim([-0.2,0.6])
            hold on
            subplot(2,2,3)
            [counts,x]=histcounts(trackArray(trackArray(:,8)==0,1));
            histogram(counts,'FaceColor',cmap(round(64*(colorInd)/maxColor),:));
            xlabel('Number of unlinked trajectories/cell')
            ylabel('Frequency')
            hold on
            % subplot(2,4,6)
            [counts,x]=histcounts(trackArray(trackArray(:,8)>0,1));
            histogram(counts,'FaceColor',cmap(round(64*(colorInd+1)/maxColor),:));
            xlabel('Number of linked trajectories/cell')
            ylabel('Frequency')
            hold on

        case 5
            figure('Name','Periodicity Analysis','NumberTitle','off');
            
            subplot(2,3,1)
            [counts,x,plotH]=KDFplotPeaks(trackArray(:,2),0.7);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            xlim([0,max(x)])
            title('Stoichiometry Distribution')
            
            subplot(2,3,2)
            [counts,x,plotH]=KDFplotPeaks(trackArray(trackArray(:,2)<(20*Isingle),2),0.7);
            set(plotH,'Color',cmap(round(64*(colorInd)/maxColor),:))
            xlabel('Stoichiometry')
            ylabel('Probability Density')
            xlim([0,20])
            title('Stoichiometry Distribution up to 20')
            
            
            subplot(2,3,3)
            PwD=pdist(trackArray(:,2));
            [counts, x]=hist(PwD,0.5:0.5:max(PwD));
            bar(x,counts)
            title('Pairwise distance distribution of stoichiometries')
            xlabel('Molecule Step')
            ylabel('Probability Density')
            %xlim([0,locs(pks==max(pks))*5])
            
            subplot(2,3,5)
            [power_spectrum_x power_spectrum_y spectrum_peaks_x spectrum_peaks_y] = FourierAndFindPeaks(x,counts,0);
            plot(power_spectrum_x, power_spectrum_y)
            xlabel('Molecule Step')
            ylabel('Power')
            xlim([0,20])
            title('Fourier Spectrum Stoichiometry up to 20')
            %xlim([0,locs(pks==max(pks))*5])
            subplot(2,3,4)
            plot(power_spectrum_x, power_spectrum_y)
            xlabel('Molecule Step')
            ylabel('Power')
            xlim([0,max(x)])
            title('Fourier Spectrum Stoichiometry')
            subplot(2,3,6)
            plot(power_spectrum_x, power_spectrum_y)
            xlabel('Molecule Step')
            ylabel('Power')
            xlim([0,10])
            title('Fourier Spectrum Stoichiometry up to 10')
        case 6
            [fitresult1,fitresult2]=CDFmobility2(trackArray(:,4));
     
    
            

        otherwise
    end
    
end
function [h,Dens, x]=KDFplotHlocal(data,bandwidthval)
colorVal=rand(1,3);
try
if exist('bandwidthval')==0
    [Dens,x] = ksdensity(data,'npoints',10000);
    h=plot(x,Dens,'color',colorVal,'LineWidth',2);
else
    [Dens,x] = ksdensity(data,'npoints',10000,'bandwidth',bandwidthval);
    h=plot(x,Dens,'color',colorVal,'LineWidth',2);
    
end
catch
    h=plot(0,0);
    x=[0,1];
    Dens=[0,1];
end
end

function AvData=binDataLocal(X,Y,numLim,num,Colour)
if nargin<3
    numLim=0;
end
if nargin<5
    Colour=rand(1,3);
end
    if nargin<4
        [counts, edges,index]=histcounts(X);
    else
        [counts, edges,index]=histcounts(X,num);
    end
for n=1:max(index)
    if sum(index==n)>0
meanY(n)=mean(Y(index==n));
%meanX(n)=mean(X(index==n));
meanX(n)=min(edges)+n*mean(edges(2:end)-edges(1:end-1));
stdY(n)=std(Y(index==n));
numY(n)=length(Y(index==n));
errorY(n)=stdY(n)/numY(n)^0.5;
    else
        meanY(n)=0;
        meanX(n)=0;
stdY(n)=0;
numY(n)=0;
errorY(n)=0;
    end
end

errorbar(meanX(numY>numLim), meanY(numY>numLim),errorY(numY>numLim),'LineWidth',2,'color',Colour);
AvData=[meanX(numY>numLim)', meanY(numY>numLim)',errorY(numY>numLim)',numY(numY>numLim)'];

end