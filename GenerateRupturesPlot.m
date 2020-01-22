%% Synthesize an SMP Trace
% using radnomly distributed chains of elastic ruptures specified by the 
% basic snow properties delta,L, and  f
clear, close all; clc;
%% Create the Saw tooth
z = 0:.001:1.5;
n = 100;
d = 1;
f = .5;
L = 10.*n;

e = linspace(0,f,n)';
F = [zeros(L/2-n/2,1);e;zeros(L/2-n/2,1)];

%% Generate Ruptures
l = linspace(0.5,1.5,25);
d = linspace(0.15,0.25,25);
f = linspace(0.075,0.15,25);
z = 0:.001:1.5;
nRuptures = 25;
lSamp = randsample(l,1);
dSamp = randsample(d,nRuptures,'true');
fSamp = randsample(f,nRuptures,'true');
zSamp = randsample(z,nRuptures,'true');
nSamp = 101;
for ii = 1:nRuptures
    R(:,ii) = [linspace(0,fSamp(ii),nSamp - 1),0]';
    Z(:,ii) = [linspace(zSamp(ii),zSamp(ii)+dSamp(ii),nSamp)]';
end
% Sort the Random Array by minimum Z coordinate
[~,sortIx] = sort(min(Z));
R = R(:,sortIx);
Z = Z(:,sortIx);
% Repeat Rupture Characteristic past L
Lix = find(max(Z)>=lSamp,1);
R(:,Lix:nRuptures) = R(:,1:(nRuptures-Lix)+1);
Z(:,Lix:nRuptures) = Z(:,1:(nRuptures-Lix)+1)+lSamp;

cumR = R(:,1);
cumZ = Z(:,1);
% Sum the Overlapping Ruptures
ii = 2;
noRuptures = 0;
while ii <= nRuptures
    % Overlapping Indecies of ii-1 iteration
    ix1 = find(Z(1,ii)-cumZ(:)<0 & Z(end,ii)-cumZ(:)>0);
    % Indices of iith iteration contained in ii-1 iteration
    ix2 = find(Z(:,ii) <= cumZ(end));
    % Indicies of iith iteration beyond ii-1 iteration
    ix3 = find(Z(:,ii) >cumZ(end));
    if noRuptures
        cumR = [cumR;R(:,ii)];
        cumZ = [cumZ;Z(:,ii)];
        noRuptures = 0;
    elseif isempty(ix1)
        tmpZ = [linspace(cumZ(end),Z(1,ii),nSamp)]';
        tmpF = zeros(length(tmpZ),1);
        cumZ = [cumZ;tmpZ];
        cumR = [cumR;tmpF];
        Z = [Z(:,1:ii-1),[linspace(min(tmpZ),max(tmpZ),nSamp)]',Z(:,ii:end)];
        R = [R(:,1:ii-1),zeros(nSamp,1),R(:,ii:end)];
        nRuptures = size(Z,2);
        noRuptures = 1; % Callback
        if ii < Lix
            Lix = Lix+1;
            disp('lix')
        end
    elseif ~noRuptures
        % Resample and sum overlapping rupture forces
        tmpF = (interp1(cumZ(ix1),cumR(ix1),Z(ix2,ii),'linear','extrap'));
        sumF = cumR(ix2) + tmpF;
        tmpZ = [linspace(cumZ(ix1(1)),cumZ(ix1(end)),length(tmpF))]';
        % Append Summed Contribution to individual ruptures
        if isempty(ix3)
            cumR = [cumR(1:(ix1(1) - 1));sumF;R(ix1(end)+1:size(R,2),ii)];
            cumZ = [cumZ(1:ix1(1)-1);tmpZ;Z(ix1(end)+1:size(Z,2),ii)];
        else
            cumR = [cumR(1:(ix1(1) - 1));sumF;R(ix3,ii)];
            cumZ = [cumZ(1:ix1(1)-1);tmpZ;Z(ix3,ii)];
        end
    end
    ii = ii +1;
end
%% Make Figures
% tmpIx = find(all(R(:,[1:Lix])==0))
% tmp = load('syntheticRuptures.mat')
% R = tmp.R; Z = tmp.Z; cumR = tmp.cumR; cumZ = tmp.cumZ;
% Lix = find(all(R(:,2:end) == R(:,1)) & Z(1,2:end) > lSamp)+1;
Lix = 15;
tmpR = [R(:,1);R(:,Lix)];
tmpZ = [Z(:,1);Z(:,Lix)];
figure();
subplot(1,3,1)
plot(tmpR,tmpZ,'k','linewidth',2)
ylim([0,2.5])
xlim([0 .2])
ylabel('Depth (mm)')
xlabel('Rupture Force (N)')
title('Element Chain')
set(gca,'fontweight','bold','fontname','serif','fontsize',11)
axis ij
subplot(1,3,2)
plot(R(:,:),Z(:,:),'k','linewidth',2)
ylim([0,2.5])
xlim([0 .2])
title('Individual Ruptures')
xlabel('Rupture Force (N)')
set(gca,'fontweight','bold','fontname','serif','fontsize',11)
axis ij
subplot(1,3,3)
plot(cumR,cumZ,'k','linewidth',2)
ylim([0,2.5])
xlim([0 1])
title('Summed Ruptures')
xlabel('Penetration Force (N)')
set(gca,'fontweight','bold','fontname','serif','fontsize',11)
axis ij
