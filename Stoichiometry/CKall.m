function filteredI=CKall(spots,bandwidth,weightVal,showOutput)

for q=1:max(spots(:,10))
    lengthSpot=size(spots(spots(:,10)==q,9),1);
    SpotI=spots(spots(:,10)==q,5);
    
    %    SpotStoich(i)=mean(SpotI(1:end))/2500;
    %if max(SpotI)>(5*Isingle) && max(SpotI)<(8*Isingle)
    try
        [XX,TX,DX,SD,DSD,XPRE]=ckfiltb2original(spots(spots(:,10)==q,5),bandwidth,weightVal);
        AllXX=[AllXX,XX'];
        %
        if showOutput==1
            figure
            %  plot(spots(spots(:,10)==q,5),'.','MarkerSize',2)
         plot(spots(spots(:,10)==q,5),'lineWidth',2)
           hold on
           scatter(1:length(XX),XX(1:end),'s','filled')
        end
    catch
        %    disp('error')
    end
    %end
end
filteredI=AllXX;

%PwD=pdist(AllXX(1:end)');
% %PwD=pdist(AllXX(AllXX<7000)');
%[counts, x]=hist(PwD,1:max(PwD));
% %maxXvalue=20;
% [power_spectrum_x power_spectrum_y spectrum_peaks_x spectrum_peaks_y] = FourierAndFindPeaks(x,counts,1,2000);
