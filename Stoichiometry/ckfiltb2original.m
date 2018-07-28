%chfilt(X,W,R).m, based on Dave Dmith's fortran routine, 03-08-2000
%window size W, power weighting factor R, 1<<R<100
%actmot/edge
%C     Chung-Kennedy filter (adapted from sed.f).
%C---------------------------------------------------------------------
%C     General formulation for linear causal+acausal filters 
%C     coupled to local variances with output Y(t) proportional 
%C     to "student-t" statistic (ANOVA). Edges identified by maxima
%C     in |Y(t)| > Yc where threshold Yc is set by user.  
%C
%C            Y(t) = (XP(t)-XM(t)/sqrt(S(t)
%C
%C     (scaled smooth derivative like student's t/sqrt(W)), with pre-step
%C     and post-step means XM,XP and variances SM,SP:
%C
%C              XM(t) = int(0,W) f(t')x(t-t')dt'
%C              XP(t) = int(0,W) f(t')x(t+t')dt'
%C              SM(t) = int(0,W) f(t')(x(t-t')-XM(t))**2 dt'
%C              SP(t) = int(0,W) f(t')(x(t+t')-XP(t))**2 dt'
%C
%C     and S(t) is average of SM(t) and SP(t).
%C          The averaging function f(t) is rectangular.
%C----------------------------------------------------------------------     
%C     Mod. of Chung and Kennedy J. Neurosci. Meth. 40, 71-86 (1991): for
%C     a step at time ts:
%C     If t < ts, variance SP is large, contaminated by step.
%C     If t > ts, variance SM is large, contaminated by step.
%C     Construct weight functions GP,GM between 0 and 1 which select the
%C     uncontaminated average:
%C
%C        GM(t) = SP**r/(SP**r + SM**r)
%C        GP(t) = SM**r/(SP**r + SM**r)      (1 << r < 100).
%C
%C     Their filtered output function is 
%C
%C            XX(t) = GM(t)*XM(t) + GP(t)*XP(t)
%C
%C     but can also define an uncontaminated variance for Y(t) output by 
%C
%C            S(t) = GM(t)*SM(t)+GP(t)*SP(t).
%C
%C     sqrt(..) is close to SM for t < ts and SP for t > ts, is the 
%C     uncontaminated variance.


function [XX,TX,DX,SD,DSD,XPRE]=ckfiltb2original(X,W,R);

%R=50;
%WP=100;
%W=10;

% Extend time series by W points at each end (reflection will do)
N=length(X);
Xnew(W+1:N+W)=X;
for I=1:W;
  Xnew(W+1-I)=X(I);
  Xnew(N+W+I)=X(N-I);
end
X=Xnew';
npts=N+2*W;%length of Xnew

wdiffx=zeros(npts,1);
datamx=[];
datx=zeros(N+W+1,W);

for n = 1:W;
   %datx(:,n) = X(n:N+W+n,1);
   %datamx = [datamx datx(:,n)];
   datamx(:,n) = X(n:N+W+n,1);
end   

wx = mean(datamx,2);
sx = std(datamx,0,2);%'0' flag normalizes by (length-1)
XP=wx(1:N);
XM=wx(W+1:N+W);
SDP=sx(1:N);
SDM=sx(W+1:N+W);

DSD=SDP-SDM;

SP=SDP.^2; %variance
SM=SDM.^2;
      
% Form switching functions, etc
RSP=SP.^R;
RSM=SM.^R;
GM=RSP./(RSP+RSM);%weighting factors
GP=RSM./(RSP+RSM);
if GM>=0&GM<=1&GP>=0&GP<=1
    S=GM.*SM+GP.*SP; %uncontaminated variance
    XX=GP.*XP+GM.*XM;
else
    S=SP;%select pre window if things gone awry with weighting fns
    XX=XP;
end

SD=sqrt(S);
SE=sqrt(S/W);
    
%YY=(XP-XM)./(sqrt(2)*SD);
TX=(XP-XM)./(sqrt(2)*SE);
DX=(XP-XM);
XPRE=XP;