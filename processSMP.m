%% Process SMP
clear; close all; clc
%% Use Paralell Processing
isParallel = 0;
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
workingDir = pwd;
% Establish the Directory Containing SMP data files
isUI = 1;
if isUI
    dataDir = uigetdir; % Opens User Interface to Select Driectory
else
    dataDir = 'input data directory here'; % ex. 'C:\SMP_data'
end
% Extract all .PNT files within dataDir
cd(dataDir)
listing = dir('*.PNT');
metaListing = dir('*.csv');
if isempty(metaListing)
    warning('Quality Control must be assessed for these data!!!')
    isQCplot = 1;
    metaData = [];
else
    isQCplot = 0;
    metaData = readtable(fullfile(dataDir,metaListing.name));
end
cd(workingDir)

% Plot SMP Profiles for QC
if isQCplot
    nfiles = size(listing,1);
    filename = cell(nfiles,1);
    SMP = cell(nfiles,1);
    for ii = 1:nfiles
        filename{ii} = [listing(ii).folder,'\',listing(ii).name];
        SMPdata = loadSMP(filename{ii});
        SMP{ii} = SMPdata;
        % Pick Snow Surface
        isPickSurface = 0;
        isDepthCorrection = 0;
        % QC Plot Force Depth Data
        h = figure(ii);
        plot(SMP{ii}.force,SMP{ii}.zF,'k'); axis ij
        set(h,'WindowStyle','docked')
        if isPickSurface
            if isDepthCorrection
                hold on; plot(SMP{ii}.force(1), SMP{ii}.zF(1), 'rx','markersize',10)
            else
                hold on; plot(SMP{ii}.force(surfIx), surfZ, 'rx','markersize',10)
            end
        end
    end
    metaData = qcSMP(listing);
    if nfiles > 1
        nametag = [listing(1).name(1:end-4),'-',listing(end).name(1:end-4)];
    else
        nametag = listing(1).name(1:end-4);
    end
    cd(dataDir)
    writetable(metaData,[nametag,'_metadata','.csv'])
    cd(workingDir)
    clear("SMP","SMPdata","filename")
end
    % Pick Snow Surface
    isPickSurface = 1;
    isDepthCorrection = 1;

%Parse Meta Data
goodSMPix = find(strcmp(metaData.Quality,'Drift') ...
    | strcmp(metaData.Quality,'Good') ...
    | strcmp(metaData.Quality,'Slight Drift'));

% Remove 'Bad' or 'Dry Run' Quality Files
files = listing(goodSMPix);
metaData = metaData(goodSMPix,:);
driftSMPix = find(strcmp(metaData.Quality,'Drift') ...
    | strcmp(metaData.Quality,'Slight Drift'));
filename = cell(length(files),1);
SMP = cell(length(files),1);
%sort Files by datetime
if ~isempty(files)
    files = nestedSortStruct(files,'date');
end
nfiles = size(files,1);
if nfiles > 1
    nametag = [files(1).name(1:end-4),'-',files(end).name(1:end-4)];
else
    nametag = files(1).name(1:end-4);
end
for ii = 1:nfiles
      filename{ii} = [files(ii).folder,'\',files(ii).name];
      SMPdata = loadSMP(filename{ii});
      SMP{ii} = SMPdata;      
      if ~isempty(metaData)
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
          surfZ = SMP{ii}.zF(surfIx);
          SMP{ii}.surface_z = surfZ;
          
          % Drift Correction
          if ismember(ii,driftSMPix)
              % Plot Drift Correction
              %             h = figure(ii*100);clf
              %             subplot(1,2,1)
              %           plot(SMP{ii}.force,SMP{ii}.zF,'k'); axis ij
              %           title('SMP Force Profile')
              %           xlabel('Force [N]')
              %           ylabel('Depth [mm]')
              %           set(gca,'fontsize',14,'fontweight','bold')
              surfIx = SMP{ii}.surface_ix - taperIx;
              zeroIx = taperIx;
              driftIx = zeroIx:surfIx;
              
              % OLS Linear Fit
              G = [ones(length(driftIx),1),SMP{ii}.zF(driftIx)'];
              d = SMP{ii}.force(driftIx);
              m = G\d; %m(1) is intercept m(2) is slope
              driftF = m(2).*SMP{ii}.zF + m(1); % Calculate Drift
              SMP{ii}.force = (SMP{ii}.force' - driftF); % Apply Drift Correction
              
              % Plot Drift Correction
              %           subplot(1,2,2)
              %           plot(SMP{ii}.force,SMP{ii}.zF,'k'); axis ij
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
              SMP{ii}.zF = SMP{ii}.dzF:SMP{ii}.dzF:(SMP{ii}.dzF*SMP{ii}.fsamp);
          end
      end
      end       
end

%% Invert SMP profile
isInvertSMP = 1; % Toggle if you wish to invert the SMP signal
if isempty(metaData)
    isInvertSMP = 0;
end
if isInvertSMP
tic
fthresh = 0.014; % Default
pfthresh = 0.1;
winsize = 10;
dz = 1; % measurements to increment by
invSMP = cell(length(files),1);
% Parallel Inversion
if isParallel
    parfor (ii = 1:length(files),nWorkers)
        invSMP{ii}=invertSMP_profile3c(SMP{ii}.zF,SMP{ii}.force,winsize,dz,pfthresh,fthresh);
    end
else
    for ii = 1:length(files)
        invSMP{ii}=invertSMP_profile3c(SMP{ii}.zF,SMP{ii}.force,winsize,dz,pfthresh,fthresh);
    end
end
    

% Export SMP Data and Inversion Results
isExport = 0;
if isExport
    cd(dataDir)
    save([nametag,'_data'],'SMP','-V7.3')
    save([nametag,'_results'],'invSMP','-V7.3')
    save([nametag,'_files'],'files','-V7.3')
    cd(workingDir)
end
toc
disp(' ')

else
    % Export SMP Data and if no Inversion Results
    isExport = 0;
    if isExport
        cd(dataDir)
        save([nametag,'_rawdata'],'SMP','-V7.3')
        save([nametag,'_files'],'files','-V7.3')
        cd(workingDir)
    end
end