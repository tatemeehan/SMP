function r=invertSMP_profile3c(z,F,winsize,dz,pfthresh,fthresh)
% HPM 06/18/08, 07/13/08 - provide fthresh again, and output Nm,N_T,Na
% NOTE: this version allows input of F,z vectors instead of pnt file
% INPUT: z = distance vector
%            F = force vector
%        winsize = size of window to process
%        dz = step size of calculation
%        pfthresh = dynamic rupture force threshold (percentage of max)
%         fthresh = static rupture force threshold
% OUTPUT: r.z = depth vector [n,1]
%         r.M = matrix with n depth values for m variables, p=2 is median,
%               p=1 and p=3 contain 95% of predictions [n,m,p]
%         r.vars = variable names [1,m]
%         r.M2 = matrix with values that do not have uncertainty
%               predictions
%         r.vars2 = variable names for M2
if nargin<5
    fthresh=0.014; % based on Lutz, pers comm, 07/08
    if nargin<4
        pfthresh=0.1; % 10% gave reasonable results, but relative results insensitive
        if nargin<3
            dz=1; % step size for calculation
        end
    end
end
%pfthresh=0.1; % dynamic rupture force threshold (noise removal)
%fthresh=0.014; % static threshold
mu=0.25; % coefficient of friction
A=pi*2.5^2;
dzF=1/250;
theta=30*pi/180; % 1/2 angle of cone
n1=floor((length(z)-winsize/dzF+1)/(dz/dzF)); % number of samples
% n1=floor((length(z)-2499)/625); % number of 10mm sections with 75% overlap
% n2=floor((length(z)-249)/250); % number of 1mm sections with 50% overlap
% preallocate
z1=zeros(n1,1); L1=zeros(n1,3); f1=zeros(n1,3); d1=zeros(n1,1); Ne1=zeros(n1,1);
Pc1=zeros(n1,3); Pc1b=zeros(n1,3); k1=zeros(n1,3); E1=zeros(n1,3); S1=zeros(n1,3);
rho1=zeros(n1,1); TI1=zeros(n1,1); FS1=zeros(n1,4); Nm1=zeros(n1,1); N_T1=zeros(n1,3);
Na1=zeros(n1,3);
% first lets process in 10 mm sections with 50% overlap
%disp('process winsize mm sections')
winsize=floor(winsize/dzF); %number of measuremnts
wstep=floor(dz/dzF); % number of measurements to increment by
%wstep=floor(winsize/4);
for n=1:n1%floor((length(z)-winsize)/wstep)
    s1=1+(n-1)*wstep; % start at 1, increment by wstep
    s2=s1+winsize-1; % process winsize mm  at a time
    r = invertSMP5(F(s1:s2),z(s1:s2),pfthresh,fthresh,mu,A,theta);
    F3=F(s1:s2); % now calculate statistical relationships
    if median(F3)>0 % make sure mean force is >0 to avoid complex values
        rho1(n)=(55.6*log(median(F3))+317.4)'; % statistical density model
    else
        rho1(n)=NaN;
    end
    TI1(n)=(1.45+5.72*std(F3)./mean(F3))'; % texture index = mean grain size/density model
    z1(n)=mean([z(s1) z(s2)]);
    L1(n,:)=r.L;
    f1(n,:)=r.fn;
    d1(n)=r.delta;
    Ne1(n)=r.Ne;
    Nm1(n)=r.Nm;
    N_T1(n,:)=r.N_T;
    Na1(n,:)=r.Na;
    Pc1(n,:)=r.Pc;
    Pc1b(n,:)=r.Pc2;
    k1(n,:)=r.k;
    E1(n,:)=r.E;
    S1(n,:)=r.sig;
    FS1(n,:)=[mean(F3) median(F3) std(F3) iqr(F3)];
end
clear r
r.z=z1(:);
M(:,:,1)=[L1(:,1) f1(:,1) Pc1(:,1) Pc1b(:,1) k1(:,1) E1(:,1) S1(:,1) N_T1(:,1) Na1(:,1) FS1(:,1)-FS1(:,3) FS1(:,2)-FS1(:,4)/2];
M(:,:,2)=[L1(:,2) f1(:,2) Pc1(:,2) Pc1b(:,2) k1(:,2) E1(:,2) S1(:,2) N_T1(:,2) Na1(:,2) FS1(:,1) FS1(:,2)];
M(:,:,3)=[L1(:,3) f1(:,3) Pc1(:,3) Pc1b(:,3) k1(:,3) E1(:,3) S1(:,3) N_T1(:,3) Na1(:,3) FS1(:,1)+FS1(:,3) FS1(:,2)+FS1(:,4)/2];
r.M=M;
r.vars={'L','f','Pc','Pc2','k','E','S','N_T','Na','mF','medF'};
r.M2=[d1(:) Ne1(:) rho1(:) TI1(:) Nm1(:)];
r.vars2={'d','Ne','rho','TI','Nm'};