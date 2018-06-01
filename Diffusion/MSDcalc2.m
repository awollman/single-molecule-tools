function [MSD,tau,weightsVal]=MSDcalc2(x,y,t)


dt=mean(t(2:end)-t(1:end-1));

%loop over time intervals
for n=1:length(t)-1
 
   MSD(n)= mean((x(1+n:end)-x(1:end-n)).^2+(y(1+n:end)-y(1:end-n)).^2);
   LengthData(n)=length((x(1+n:end)-x(1:end-n)).^2+(y(1+n:end)-y(1:end-n)).^2);
    tau(n)=n*dt;
end
weightsVal=LengthData/max(LengthData);
end