clear; close all; clc
isWindows = 1;
isLinux = 0;
isParallel = 1;
if isParallel
    % Wake Parallel Computing
    if isempty(gcp('nocreate'))     
        nWorkers = 4;
        p = parpool(nWorkers);
    else
        nWorkers = 4;
    end
end
%% Read SMP Profile
Locs = 1:10; % Number of SMP Sites
Coords = cell(length(Locs),1);
for ff = Locs
    workingDir = pwd;
    % 25Jan2018 data
    if ff == 1
        %       25CRWYSMPS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/25Jan2018/SMPS/25CRWYSMPS';
        end
        if isWindows
%             dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\25Jan2018\SMPS\25CRWYSMPS';
            dataDir = 'D:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\25Jan2018\SMPS\25CRWYSMPS';
        end        
        tag = '25CRWYSMPS';
        disp(['Processing ',tag]);
        metaData = readtable(fullfile(dataDir,'25CRWYSMPSdataQC.xlsx'));
    end
    if ff == 2
        % 25CTWYSMPS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/25Jan2018/SMPS/25CTWYSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\25Jan2018\SMPS\25CTWYSMPS';
        end        
        tag = '25CTWYSMPS';
        disp(['Processing ',tag]);        
        metaData = readtable(fullfile(dataDir,'25CTWYSMPSdataQC.xlsx'));
    end
    % 26Jan2018 data
    if ff == 3
        % 26MobLPSMPS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/26Jan2018/SMPS/26MobLPSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\26Jan2018\SMPS\26MobLPSMPS';
        end          
        tag = '26MobLPSMPS';
        disp(['Processing ',tag]);        
        metaData = readtable(fullfile(dataDir,'26MobLPSMPSdataQC.xlsx'));
    end
    
    % 27Jan2018 data
    if ff == 4
        % 27SWVSSPMS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/27Jan2018/SMPS/27SWVSSPMS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\27Jan2018\SMPS\27SWVSSPMS';
        end         
        tag = '27SWVSSPMS';
        disp(['Processing ',tag]);        
        metaData = readtable(fullfile(dataDir,'27SWVSSMPSdataQC.xlsx'));
    end
    % 27SEVSSMPS
    if ff == 5
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/27Jan2018/SMPS/27SEVSSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\27Jan2018\SMPS\27SEVSSMPS';
        end           
        tag = '27SEVSSMPS';
        disp(['Processing ',tag]);        
        metaData = readtable(fullfile(dataDir,'27SEVSSPMSdataQC.xlsx'));
    end
    % 28Jan2018 data
    if ff == 6
        % 28CMobLPSMPS
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/28Jan2018/SMPS/28CMobLPSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\28Jan2018\SMPS\28CMobLPSMPS';
        end        
        tag = '28CMobLPSMPS';
        disp(['Processing ',tag]);        
        metaData = readtable(fullfile(dataDir,'28CMobLPSMPSdataQC.xlsx'));
    end
    
    % 29Jan2018 data
    if ff == 7
        % 29SWVSSPM data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/29Jan2018/SMPS/29SWVSSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\29Jan2018\SMPS\29SWVSSMPS';
        end        
        tag = '29SWVSSMPS';
        disp(['Processing ',tag]);        
        metaData = readtable(fullfile(dataDir,'29SWVSSMPSdataQC.xlsx'));
    end
    
    % 30Jan2018 data
    if ff == 8
        % 30NMobLPSPMS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/30Jan2018/SMPS/30NMobLPSPMS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\30Jan2018\SMPS\30NMobLPSPMS';
        end        
        tag = '30NMobLPSPMS';
        disp(['Processing ',tag]);        
        metaData = readtable(fullfile(dataDir,'30NMobLPSMPSdataQC.xlsx'));
    end
    
    if ff == 9
        % 30STXYMTVRCDSMPS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/30Jan2018/SMPS/30STXYMTVRCDSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\30Jan2018\SMPS\30STXYMTVRCDSMPS';
        end         
        tag = '30STXYMTVRCDSMPS';
        disp(['Processing ',tag]);
        metaData = readtable(fullfile(dataDir,'30STXYMTVRCDSMPSdataQC.xlsx'));
    end
    
    % 31Jan2018 data
    if ff == 10
        % 31STXYVSMTVRCDSMPS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/31Jan2018/SMPS/31STXYVSMTVRCDSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\31Jan2018\SMPS\31STXYVSMTVRCDSMPS';
        end        
        tag = '31STXYVSMTVRCDSMPS';
        disp(['Processing ',tag]);
        metaData = readtable(fullfile(dataDir,'31STXYVSMTVRCDSMPSdataQC.xlsx'));
    end

% Plot SMP Profiles for QC
isQCplot = 0;
% Pick Snow Surface
isPickSurface = 1;
isDepthCorrection = 1;

%Parse Meta Data
% groupSMP = findgroups(metaData.Quality); % Group Data by Quality
% Find 'Drift','Good','Slight Drift'; 1 == 'Bad'; 3 == 'Dry Run';
goodSMPix = find(strcmp(metaData.Quality,'Drift') ...
    | strcmp(metaData.Quality,'Good') ...
    | strcmp(metaData.Quality,'Slight Drift'));
% goodSMPix = find((groupSMP ~= 1) & (groupSMP ~= 3));
% groupSMP = groupSMP(goodSMPix);
% driftSMPix = find((groupSMP ~= 4));

% Get All SMP Files
files = dir(fullfile(dataDir,'*.pnt'));
headers = dir(fullfile(dataDir,'*.txt'));
configIx = find(strcmp({headers.name},'Config.txt'));
% Remove 'Bad' or 'Dry Run' Quality Files
headers(configIx) = [];
files = files(goodSMPix);
headers = headers(goodSMPix);
metaData = metaData(goodSMPix,:);
driftSMPix = find(strcmp(metaData.Quality,'Drift') ...
    | strcmp(metaData.Quality,'Slight Drift'));
filename = cell(length(files),1);
SMP = cell(length(files),1);
%sort Files by datetime
if ~isempty(files)
    files = nestedSortStruct(files,'date');
end
for ii = 1:length(files)
      filename{ii} = [files(ii).folder,'/',files(ii).name];
      SMPdata = readSMP(filename{ii});
      headerOpts = detectImportOptions(fullfile(dataDir,headers(ii).name));
      headerOpts.DataLine = 17;
      header = readtable(fullfile(dataDir,headers(ii).name),headerOpts);
      LongIx = find(strcmp(header.Var1,'Longitude'));
      LatIx = find(strcmp(header.Var1,'Latitude'));
      SMPdata.depth_F=SMPdata.dzF:SMPdata.dzF:(SMPdata.dzF*length(SMPdata.force));
      SMPdata.depth_T=SMPdata.dzT:SMPdata.dzT:(SMPdata.dzT*length(SMPdata.temp));
      SMPdata.longitude = str2num(cell2mat(header.Var2(LongIx)));
      SMPdata.latitude = str2num(cell2mat(header.Var2(LatIx)));
      SMP{ii} = SMPdata;
      Coords{ff}(ii,:) = [SMPdata.longitude,SMPdata.latitude];
      
      % First Break Pick Linear Drift Correction
      if isPickSurface
          filtWin = 250;
          trace = SMP{ii}.force; % Query Signal
          trace = movmean(trace,filtWin/4);% Smooth Trace
          taperIx = 2*filtWin;
          vartrace = movvar(trace,filtWin);% Estimate Signal Variance
          airvar = mean(vartrace(taperIx:(3*taperIx)));
          surfIx = find(vartrace(taperIx:end) > 49.*airvar, 1)+taperIx;% 7 sigma Threshold
          SMP{ii}.surface_ix = surfIx;
          surfZ = SMP{ii}.depth_F(surfIx);
          SMP{ii}.surface_z = surfZ;
          
          % Drift Correction
          if ismember(ii,driftSMPix)
              % Plot Drift Correction
              %             h = figure(ii*100);clf
              %             subplot(1,2,1)
              %           plot(SMP{ii}.force,SMP{ii}.depth_F,'k'); axis ij
              %           title('SMP Force Profile')
              %           xlabel('Force [N]')
              %           ylabel('Depth [mm]')
              %           set(gca,'fontsize',14,'fontweight','bold')
              surfIx = SMP{ii}.surface_ix - taperIx;
              zeroIx = taperIx;
              driftIx = zeroIx:surfIx;
              
              % OLS Linear Fit
              G = [ones(length(driftIx),1),SMP{ii}.depth_F(driftIx)'];
              d = SMP{ii}.force(driftIx);
              m = G\d; %m(1) is intercept m(2) is slope
              driftF = m(2).*SMP{ii}.depth_F + m(1); % Calculate Drift
              SMP{ii}.force = (SMP{ii}.force' - driftF); % Apply Drift Correction
              
              % Plot Drift Correction
              %           subplot(1,2,2)
              %           plot(SMP{ii}.force,SMP{ii}.depth_F,'k'); axis ij
              %           hold on; plot(SMP{ii}.force(surfIx), surfZ, 'rx','markersize',24,'linewidth',2)
              %           title('Drift Corrected Profile')
              %           xlabel('Force [N]')
              %           ylabel('Depth [mm]')
              %           set(gca,'fontsize',14,'fontweight','bold')
              %           legend('Force','Snow Surface')
          end
          % Correct Penetration Depth
          if isDepthCorrection
              SMP{ii}.force(1:SMP{ii}.surface_ix-1) = [];
              SMP{ii}.fsamp = length(SMP{ii}.force);
              SMP{ii}.depth_F = SMP{ii}.dzF:SMP{ii}.dzF:(SMP{ii}.dzF*SMP{ii}.fsamp);
          end
      end
          
      % QC Plot Force Depth Data
      if isQCplot
          h = figure(ii);
          plot(SMP{ii}.force,SMP{ii}.depth_F,'k'); axis ij
          if isPickSurface
              if isDepthCorrection
              hold on; plot(SMP{ii}.force(1), SMP{ii}.depth_F(1), 'rx','markersize',10)
              else
              hold on; plot(SMP{ii}.force(surfIx), surfZ, 'rx','markersize',10)
              end
              set(h,'WindowStyle','docked')
          end
      end
     
end
% Get Average Location of Acquisition for Waypoint Positions 
try
    Waypoints(ff).location = tag;
    Waypoints(ff).longitude = mean(Coords{ff}(:,1));
    Waypoints(ff).latitude = mean(Coords{ff}(:,2));
    Waypoints(ff).coords =  Coords{ff};
    Waypoints(ff).files = {files.name};
catch
    Waypoints(ff).longitude = [];
    Waypoints(ff).latitude = [];
end
%% Invert SMP profile
isInvert = 1;
if isInvert
fthresh = 0.0195; % Default
pfthresh = 0.1;
winsize = 10;
dz = 1; % [mm] measurements to increment by
invSMP = cell(length(files),1);
% Parallel Inversion
if isParallel
    parfor (ii = 1:length(files),nWorkers)
        invSMP{ii}=invertSMP_profile3c(SMP{ii}.depth_F,SMP{ii}.force,winsize,dz,pfthresh,fthresh);
    end
else
    for ii = 1:length(files)
        invSMP{ii}=invertSMP_profile3c(SMP{ii}.depth_F,SMP{ii}.force,winsize,dz,pfthresh,fthresh);
    end
end
    

% Export SMP Data and Inversion Results
isExport = 1;
if isExport
    cd(dataDir)
    save([tag,'data'],'SMP','-V7.3')
    % save([tag,'results'],'invSMP','-V7.3')
    % save([tag,'files'],'files','-V7.3')
    cd(workingDir)
end
toc
disp(' ')
% Clear Memory for next Loop
clear
else
    % Export SMP Data and if no Inversion Results
    isExport = 0;
    if isExport
        cd(dataDir)
        save([tag,'data'],'SMP','-V7.3')
        % save([tag,'results'],'invSMP','-V7.3')
        % save([tag,'files'],'files','-V7.3')
        cd(workingDir)
    end
end
% Also Export the Waypoints
if isExport
    cd('E:\CRREL_SnowCompaction\W_YELLOWSTONE\GPS')
    save('SMPwaypoints','Waypoints','-V7.3')
    cd(workingDir)
end
end

%% Estimate Noise Level Force Threshold
% Estimate fthresh from data
isCalFthresh = 0;
if isCalFthresh
tic
c = 1;
for jj = [1,3,5,10,25,50,125,250]
    c
    noiseRange = [];
for ii = 1:length(files)
    
surfIx = SMP{ii}.surface_ix - taperIx;
zeroIx = taperIx;
noiseIx = zeroIx:surfIx;
noiseMax = movmax(SMP{ii}.force(noiseIx),jj);
noiseMin = movmin(SMP{ii}.force(noiseIx),jj);
noiseMax = noiseMax(:);
noiseMin = noiseMin(:);
noiseRange = [noiseRange;(noiseMax-noiseMin)];

end

% figure();boxplot(noise,'Labels',{'1','3','5','10','25','50','125','250'})
% set(gca,'xtick',[1,3,5,10,25,50,125,250]);
%% Invert SMP Profile
% r = invertSMP5(SMP.force,SMP.depth_F,0.01);%,fthresh,mu,A,theta)
% r=invertSMP_profile3c(SMPdata.depth_F,SMPdata.force,winsize,SMPdata.dzF,pfthresh,fthresh)
% Estimate Rupture Force Threshold .95 quantile of Instrument Noise Level 
fthresh = quantile(noiseRange,0.5);
% fthresh = 0.014; % Default
pfthresh = 0.1;
winsize = 1;
invSMP = cell(length(files),1);
% Parallel Inversion
nWorkers = 12;
% parfor (ii = 1:length(files))
for kk = 18
invSMP{c}=invertSMP_profile3c(SMP{kk}.depth_F,SMP{kk}.force,winsize,SMP{kk}.dzF,pfthresh,fthresh);

h=figure(jj);plot(SMP{kk}.force,SMP{kk}.depth_F,'k')
axis ij
hold on;
plot(SMP{kk}.force(find(isnan(invSMP{c}.M2(:,1)))),SMP{kk}.depth_F(find(isnan(invSMP{c}.M2(:,1)))),'r.')
legend('Force','No Ruptures')
title(['Force Threshold: ',num2str(fthresh)])
xlabel('Force [N]')
ylabel('Uncorrected Depth [mm]')
set(gca,'fontweight','bold','fontsize',14)
set(h,'WindowStyle','docked')
end
noise(:,c) = noiseRange;
    c = c+1;
end
h2=figure(999);boxplot(noise,'Labels',{'1','3','5','10','25','50','125','250'})
set(gca,'xtick',[1,3,5,10,25,50,125,250]);
title('Rupture Force Threshold')
xlabel('Window Length')
ylabel('Noise Floor [N]')
set(gca,'fontweight','bold','fontsize',14)
set(h2,'WindowStyle','docked')
toc
end