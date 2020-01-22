%% SMP Comparison with Shear Vane
% Tate Meehan 6/5/18
clear; close all; clc;
isWindows = 1;
isLinux = 0;
%% Load SMP Data and Inversion Result
locs = [1,2,3,4,5,6,7,8,9,10,11];
 workingDir = pwd;
 medM = cell(length(locs),1);
 nn = 0;
%  shearVane = [];
shearVane = cell(length(locs),1);
SVtag = cell(1,1);
q = 0.5; % Test Quantile for Correlations
depth = 50;
for ff = locs
    % 25Jan2018 data
    if ff == 1
        %       25CRWYSMPS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/25Jan2018/SMPS/25CRWYSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\25Jan2018\SMPS\25CRWYSMPS';
        end
        
        tag = '25CRWYSMPS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        SVtag{1} = [SVtag{1},{'25CRWY'}];
        % Shear Vane Depth
%         depth = 50;
        % Runway Ruts
        nn = nn + 1;
        % Rut files
        rutVane = [73,93,80,51,52,88,83,98,92];
        rutFiles = [45:64];
        badIx = find(rutFiles == 51);%rmv 51
        rutFiles(badIx) = [];
        rutIx = rutFiles;
        shearVane{nn} = rutVane;
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = rutIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
%         medM{nn}( find(medM{nn}), :, : ) = [];
%         medM{nn}( ~any(medM{nn},2), : ) = [];
    end
    if ff == 2
        % 25CTWYSMPS data
        tic
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\25Jan2018\SMPS\25CTWYSMPS';
        end
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/25Jan2018/SMPS/25CTWYSMPS';
        end
        tag = '25CTWYSMPS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        SVtag{1} = [SVtag{1},{'25CTWY'}];
        
        % Shear Vane Depth
%         depth = 50;
        % Runway Ruts
        nn = nn + 1;
        % Rut files
        rutVane = [72,48,6,130,55];
        rutFiles = [16:30];
        %         badIx = find(rutFiles == 51);%rmv 51
        %         rutFiles(badIx) = [];
        rutIx = rutFiles;
        shearVane{nn} = rutVane;
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = rutIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
        %         medM{nn}( ~any(medM{nn},2), : ) = [];
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
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
                
        % Mobility Loop Rut 1
        SVtag{1} = [SVtag{1},{'26MobLP-1'}];
        nn = nn + 1;
        % Rut files
        rutVane = [55,88,38,52,56,48];
        rutFiles = [254,255,256,257,258,259,264,265,267,268];
        
%         depth = 50; % [mm]
        
        % Loop over rut Files
        rutIx = [];
        shearVane{nn} = rutVane;
        active = 1;
        for ii = 1:length(files)
            tmp = strcmp(files(ii).name(:,6:8),num2str(rutFiles(active)));
            if tmp; rutIx = [rutIx,ii]; active = active + 1; end
            if active > length(rutFiles); active = length(rutFiles); end
        end
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = rutIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
        
%         medM{nn}( ~any(medM{nn},2), : ) = [];
    end
    
    if ff == 4
        % 26MobLPSMPS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/26Jan2018/SMPS/26MobLPSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\26Jan2018\SMPS\26MobLPSMPS';
        end         
        tag = '26MobLPSMPS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        
        
        % Mobility Loop Rut 2
        SVtag{1} = [SVtag{1},{'26MobLP-2'}];
        nn = nn + 1;
        % Rut files
        rutVane = [32,72,72,56,42];
        rutFiles = [242,243,244,245,246,247,248];
        
%         depth = 50; % [mm]
        
        % Loop over rut Files
        rutIx = [];
        shearVane{nn} = rutVane;
        active = 1;
        for ii = 1:length(files)
            tmp = strcmp(files(ii).name(:,6:8),num2str(rutFiles(active)));
            if tmp; rutIx = [rutIx,ii]; active = active + 1; end
            if active > length(rutFiles); active = length(rutFiles); end
        end
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = rutIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
%         medM{nn}( ~any(medM{nn},2), : ) = [];
    end
    if ff == 5
        % 26MobLPSMPS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/26Jan2018/SMPS/26MobLPSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\26Jan2018\SMPS\26MobLPSMPS';
        end         
        tag = '26MobLPSMPS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        
        % Mobility Loop Rut 3
        SVtag{1} = [SVtag{1},{'26MobLP-3'}];
        nn = nn + 1;
        % Rut files
        rutVane = [40,22,12];
        rutFiles = [232,233,234,238,239,241];
        
%         depth = 50; % [mm]
        
        % Loop over rut Files
        rutIx = [];
        shearVane{nn} = rutVane;
        active = 1;
        for ii = 1:length(files)
            tmp = strcmp(files(ii).name(:,6:8),num2str(rutFiles(active)));
            if tmp; rutIx = [rutIx,ii]; active = active + 1; end
            if active > length(rutFiles); active = length(rutFiles); end
        end
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = rutIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
        
%         medM{nn}( ~any(medM{nn},2), : ) = [];
    end
    if ff == 6
        % 26MobLPSMPS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/26Jan2018/SMPS/26MobLPSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\26Jan2018\SMPS\26MobLPSMPS';
        end         
        tag = '26MobLPSMPS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        % Mobility Loop Rut 4
        SVtag{1} = [SVtag{1},{'26MobLP-4'}];
        nn = nn + 1;
        % Rut files
        rutVane = [74,46,88,42];
        rutFiles = [222,223,224,225,226,227];
        
%         depth = 50; % [mm]
        
        % Loop over rut Files
        rutIx = [];
        shearVane{nn} = rutVane;
        active = 1;
        for ii = 1:length(files)
            tmp = strcmp(files(ii).name(:,6:8),num2str(rutFiles(active)));
            if tmp; rutIx = [rutIx,ii]; active = active + 1; end
            if active > length(rutFiles); active = length(rutFiles); end
        end
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = rutIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
        
%         medM{nn}( ~any(medM{nn},2), : ) = [];
    end
    % No Shear Vane Data Here!
    %     % 27Jan2018 data
    %     if ff == 4
    %         % 27SWVSSPMS data
    %         tic
    %         dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/27Jan2018/SMPS/27SWVSSPMS';
    %         tag = '27SWVSSPMS';
    %         load(fullfile(dataDir,[tag,'data']));
    %         load(fullfile(dataDir,[tag,'results']));
    %         load(fullfile(dataDir,[tag,'files']))
    %     end
    %     % 27SEVSSMPS
    %     if ff == 5
    %         tic
    %         dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/27Jan2018/SMPS/27SEVSSMPS';
    %         tag = '27SEVSSMPS';
    %         load(fullfile(dataDir,[tag,'data']));
    %         load(fullfile(dataDir,[tag,'results']));
    %         load(fullfile(dataDir,[tag,'files']))
    %     end
    %     % 28Jan2018 data
    %     if ff == 6
    %         % 28CMobLPSMPS
    %         tic
    %         dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/28Jan2018/SMPS/28CMobLPSMPS';
    %         tag = '28CMobLPSMPS';
    %         load(fullfile(dataDir,[tag,'data']));
    %         load(fullfile(dataDir,[tag,'results']));
    %         load(fullfile(dataDir,[tag,'files']))
    %     end
    
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
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        
        % Razor CD Ruts
        SVtag{1} = [SVtag{1},{'29SWVS-RZR'}];
        nn = nn + 1;
        % Rut files
        rutVane = [58,24,62,20,30,78,46,62,53,130];
        rutFiles = [404,405,416,417]; %rmv 404
        
%         depth = 50; % [mm]
        
        % Loop over rut Files
        rutIx = [];
        shearVane{nn} = rutVane;
        active = 1;
        for ii = 1:length(files)
            tmp = strcmp(files(ii).name(:,6:8),num2str(rutFiles(active)));
            if tmp; rutIx = [rutIx,ii]; active = active + 1; end
            if active > length(rutFiles); active = length(rutFiles); end
        end
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = rutIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
        
%         medM{nn}( ~any(medM{nn},2), : ) = [];
    end
    if ff == 8
        % 29SWVSSPM data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/29Jan2018/SMPS/29SWVSSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\29Jan2018\SMPS\29SWVSSMPS';
        end        
        tag = '29SWVSSMPS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        
        % LVSR CD Ruts
        SVtag{1} = [SVtag{1},{'29SWVS-LSVR'}];
        nn = nn + 1;
        % Rut files
        rutVane = [33,14,33,12,130,14,34,76,1,30,22,26];
        rutFiles = [440,448,449,450];
        
%         depth = 50; % [mm]
        
        % Loop over rut Files
        rutIx = [];
        shearVane{nn} = rutVane;
        active = 1;
        for ii = 1:length(files)
            tmp = strcmp(files(ii).name(:,6:8),num2str(rutFiles(active)));
            if tmp; rutIx = [rutIx,ii]; active = active + 1; end
            if active > length(rutFiles); active = length(rutFiles); end
        end
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = rutIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum        
%         medM{nn}( ~any(medM{nn},2), : ) = [];
        
    end
    
    % 30Jan2018 data
    if ff == 9
        % 30NMobLPSPMS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/30Jan2018/SMPS/30NMobLPSPMS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\30Jan2018\SMPS\30NMobLPSPMS';
        end         
        tag = '30NMobLPSPMS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        nn = nn+1;
        % Rut files
        rutVane = [62,85,130,122,24,62,130,82,81,90];
        rutFiles = [460,461,467,468,469,486];
        % Belly Files
        bellyFiles = [462,463,464,465,466];
        bellyVane = [66,89,82,70,50,78,56,24,56,20]; % reading kPa
        %         shearVane = [shearVane;bellyVane'];
        
%         depth = 50; % [mm]
        
        % Loop over rut Files
        SVtag{1} = [SVtag{1},{'30NMobLP-R'}];
        rutIx = [];
        shearVane{nn} = rutVane;
        active = 1;
        for ii = 1:length(files)
            tmp = strcmp(files(ii).name(:,6:8),num2str(rutFiles(active)));
            if tmp; rutIx = [rutIx,ii]; active = active + 1; end
            if active > length(rutFiles); active = length(rutFiles); end
        end
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = rutIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
        
%         medM{nn}( ~any(medM{nn},2), : ) = [];
    end
    
    if ff == 10
        % 30NMobLPSPMS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/30Jan2018/SMPS/30NMobLPSPMS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\30Jan2018\SMPS\30NMobLPSPMS';
        end        
        tag = '30NMobLPSPMS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        % Loop over belly Files
        SVtag{1} = [SVtag{1},{'30NMobLP-B'}];
        bellyIx = [];
        nn = nn+1;
        shearVane{nn} = bellyVane;
        active = 1;
        for ii = 1:length(files)
            tmp = strcmp(files(ii).name(:,6:8),num2str(bellyFiles(active)));
            if tmp; bellyIx = [bellyIx,ii]; active = active + 1; end
            if active > length(bellyFiles); active = length(bellyFiles); end
        end
        
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = bellyIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
        
%         medM{nn}( ~any(medM{nn},2), : ) = [];
    end
    
    % No Shear Vane Here!
    %     if ff == 9
    %         % 30STXYMTVRCDSMPS data
    %         tic
    %         dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/30Jan2018/SMPS/30STXYMTVRCDSMPS';
    %         tag = '30STXYMTVRCDSMPS';
    %         load(fullfile(dataDir,[tag,'data']));
    %         load(fullfile(dataDir,[tag,'results']));
    %         load(fullfile(dataDir,[tag,'files']))
    %     end
    
    % 31Jan2018 data
    if ff == 11
        % 31STXYVSMTVRCDSMPS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/31Jan2018/SMPS/31STXYVSMTVRCDSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\31Jan2018\SMPS\31STXYVSMTVRCDSMPS';
        end        
        tag = '31STXYVSMTVRCDSMPS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        
        % Shear Vane Data ~ MTVR Belly
        SVtag{1} = [SVtag{1},{'31SRWYVS-B'}];
        nn = nn+1;
        bellyVane = [17,8,6,24,12,16,8,0,14,10,5,3,6,14]; % reading kPa
        %         shearVane = [shearVane;bellyVane'];
        shearVane{nn} = bellyVane;
%         depth = 50; % [mm]
        
        % SMP Belly 521 - 525
        % Loop over files to check file names; inconvienient ikr
        bellyFiles = [521,522,523,524,525]; % manual entry from Susans Note
        bellyIx = [];
        
        active = 1;
        for ii = 1:length(files)
            tmp = strcmp(files(ii).name(:,6:8),num2str(bellyFiles(active)));
            if tmp; bellyIx = [bellyIx,ii]; active = active + 1; end
            if active > length(bellyFiles); active = length(bellyFiles); end
        end
        
        % Calculate Median SMP Value for Shear Vane Depth
        for ii = bellyIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);            
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
%                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
%                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
%                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
%                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);                
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
%             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.34);
%             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx2),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
%         medM{nn}( ~any(medM{nn},2), : ) = [];
    end
end
 
 % Concatenate Belly data
%  Mbelly = cell2mat(Mbelly);
%% Median Statistic of Shear Vane and SMP for Correlation
boxShearVane = [];
boxShearVaneIx = [];
boxSMP = [];
boxSMPIx = [];
for ii = 1:size(medM,1)
    medSMP(ii,:) = nanmedian(medM{ii}(:,:,1),1);
    SMP5(ii,:) = nanmedian(medM{ii}(:,:,2),1);
    SMP95(ii,:) = nanmedian(medM{ii}(:,:,3),1);
    medShearVane(ii,1) = median(shearVane{ii});
    SV5(ii,1) = quantile(shearVane{ii},0.05);
    SV95(ii,1) = quantile(shearVane{ii},0.95);
    boxShearVane = [boxShearVane;shearVane{ii}(:)];
    boxShearVaneIx = [boxShearVaneIx;ii*ones(length(shearVane{ii}),1)];
    boxSMP = [boxSMP;medM{ii}];
    boxSMPIx = [boxSMPIx;ii*ones(size(medM{ii},1),1)];
end
% Compute Correlation Coefficient
for ii = 1:size(medSMP,2)
    [R,PP,RL,RU] = corrcoef(medSMP(:,ii),medShearVane);
    CC(1,ii) = R(1,2);
    CC(2,ii) = PP(1,2);
    CC(3,ii) = RL(1,2);
    CC(4,ii) = RU(1,2);
end
%% Visualization
% Box Plots
% Shear Vane
% SVtag = {['25CRWY'],['25CTWY'],['26MobLP-1'],['26MobLP-2'],['26MobLP-3'],['26MobLP-4'],['29SWVS-RZR']...
%     ['29SWV-LSVR'],['30NMobLP-R'],['30NMobLP-B'],['31STXYVS-B']};

% figure();boxplot(boxShearVane,boxShearVaneIx)
% 
% figure();
% for ii = 1:size(boxSMP,2)
%     subplot(4,5,ii)
%     boxplot(boxSMP(:,ii),boxSMPIx)
% end

% Correlation Plot
figure();
titleVar = cat(2,[invSMP{1}.vars],[invSMP{1}.vars2],'force');
titleVar = cat(2,[invSMP{1}.vars],[invSMP{1}.vars2],'F');
titleVar{8} = 'N_{T}';titleVar{12} = '\delta';titleVar{14} = '\rho';
titleVar{10} = 'F_{m}'; titleVar{15} = 'T_{I}';titleVar{11} = 'F_{med}';
titlevar{9} = 'N_{a}';titleVar{3} = 'P_{c1}';titleVar{4} = 'P_{c2}';
titleVar{13} = 'N_e';titleVar{9} = 'N_a';titleVar{16} = 'N_{m}';
titleStr = {['Characteristic Element Length'],['Rupture Force per mm'],['Probability of Contact'],['Probability of Contact'],['Stiffness'],...
    ['Elastic Modulus'],['Compressive Strength'], ['Number of Ruptures per Characteristic Length'],['Number of Available Elements'],['Mean Rupture Force per mm'],...
    ['Median Rupture Force per mm'],['Deflection'],['Number of Engaged Elements'],['Density'],['Textural Index'],['Number of Measured Ruptures per mm'],['Penetration Force']};
subarray = [1,2,12,3,4,13,9,17,10,11,14,15,8,16,5,6,7];

for ii = 1:size(boxSMP,2)
    subplot(5,4,ii)
%     plot(medSMP(:,ii),medShearVane,'ok')
    plot(medSMP(:,subarray(ii)),medShearVane,'ok')

    hold on; lsline
    set(gca,'fontweight','bold')
%     title(titleVar{ii})
    title(titleVar{subarray(ii)})
%     text(0.675,0.875,['R = ', num2str(CC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.675,0.675,['p = ', num2str(CC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.675,0.875,['R = ', num2str(CC(1,subarray(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.675,0.675,['p = ', num2str(CC(2,subarray(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    
    if ii == size(boxSMP,2)
        subplot(5,4,[18,19,20])
        boxplot(boxShearVane,boxShearVaneIx,'Labels',SVtag)
        xtickangle(-10)
        ylabel('Shear Vane (kPa)','rotation',270,'units','normalized','position',[1.0625,.5,0])
        xlabel('Location','units','normalized','position',[.5,-.25,0])
        set(gca,'YAxisLocation','right','fontweight','bold')
%         suplabel('SMP','x')
%         suplabel('Shear Vane','y')
    end
end

% Significant Correlations
sigIx = find(CC(2,:) <0.1);
if length(sigIx)<4;
    subcols = length(sigIx);
else
    subcols  = 4;
end
subrows = ceil(length(sigIx)/subcols);

figure();
for ii = 1:length(sigIx)
    subplot(subrows + 1,subcols,ii)
    plot(medSMP(:,sigIx(ii)),medShearVane,'ok'); hold on;
    lsline;
%     plot(SMP5(:,sigIx(ii)),medShearVane,'or','markerfacecolor','r');
%     plot(SMP95(:,sigIx(ii)),medShearVane,'or','markerfacecolor','r');
%     errorbar(medSMP(:,sigIx(ii)),medShearVane,(medShearVane-SV5),-(medShearVane-SV95),(medSMP(:,sigIx(ii))-SMP5(:,sigIx(ii))),-(medSMP(:,sigIx(ii))-SMP95(:,sigIx(ii))),'ok')
    errorbar(medSMP(:,sigIx(ii)),medShearVane,(medSMP(:,sigIx(ii))-SMP5(:,sigIx(ii))),-(medSMP(:,sigIx(ii))-SMP95(:,sigIx(ii))),'horizontal','ok')

    text(0.775,0.925,['R = ', num2str(CC(1,sigIx(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.775,0.85,['p = ', num2str(CC(2,sigIx(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    title(titleVar{sigIx(ii)})
    set(gca,'fontweight','bold')
    ylim([0,135])
    ylabel('Shear Vane (kPa)')
    xlabel(titleStr{sigIx(ii)})
    
%     axis square
    % add a box plot
        if ii == length(sigIx)
        subplot(subrows + 1,subcols,[ii+1:ii+subcols])
        boxplot(boxShearVane,boxShearVaneIx,'Labels',SVtag)
        xtickangle(-10)
        ylabel('Shear Vane (kPa)','rotation',270,'units','normalized','position',[1.0625,.5,0])
        xlabel('Location','units','normalized','position',[.5,-.1,0])
        set(gca,'YAxisLocation','right','fontweight','bold')
%         suplabel('SMP','x')
%         suplabel('Shear Vane','y')
        end
end
    
%% Correlation Plot
figure();
titleVar = cat(2,[invSMP{1}.vars],[invSMP{1}.vars2],'F');
titleVar{8} = 'N_{T}';titleVar{12} = '\delta';titleVar{14} = '\rho';
titleVar{10} = 'F_{m}'; titleVar{15} = 'T_{I}';titleVar{11} = 'F_{med}';
titlevar{9} = 'N_{a}';titleVar{3} = 'P_{c1}';titleVar{4} = 'P_{c2}';
titleVar{13} = 'N_e';titleVar{9} = 'N_a';titleVar{16} = 'N_{m}';
titleStr = {['Characteristic Element Length'],['Rupture Force per mm'],['Probability of Contact'],['Probability of Contact'],['Stiffness'],...
    ['Elastic Modulus'],['Compressive Strength'], ['Number of Ruptures per Characteristic Length'],['Number of Available Elements'],['Mean Rupture Force per mm'],...
    ['Median Rupture Force per mm'],['Deflection'],['Number of Engaged Elements'],['Density'],['Textural Index'],['Number of Measured Ruptures per mm'],['Penetration Force']};
% subarray = [1,2,3,4,5,6,7,8,9,10,11,13,14,15,17,18,19];
% plotarray = [1,2,3,4,5,6,7,8,9,10,11,13,14,15,17,18,19];
plotarray = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17];
subarray = [1,2,12,3,4,13,9,17,10,11,14,15,8,16,5,6,7];
sortMedSMP = medSMP(:,subarray);
for ii = 1:size(boxSMP,2)
    subplot(3,6,plotarray(ii))
%         plot(medSMP(:,ii),medramm,'ok')

    plot(medSMP(:,subarray(ii)),medShearVane,'ok')
%     plot(SMP95(:,subarray(ii)),medramm,'ok')
%     plot(normMedSMP(:,subarray(ii)),medramm,'ok')

    hold on; lsline
    set(gca,'fontweight','bold')
%     title(titleVar{ii})
    title(titleVar{subarray(ii)})
%     text(0.675,0.875,['R = ', num2str(CC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.675,0.675,['p = ', num2str(CC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.075,0.775,['R = ', num2str(XC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.075,0.575,['p = ', num2str(XC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
    if CC(1,subarray(ii)) < 0
    text(0.525,0.875,['R = ', num2str(CC(1,subarray(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.525,0.775,['p = ', num2str(CC(2,subarray(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    else 
        text(0.525,0.275,['R = ', num2str(CC(1,subarray(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.525,0.175,['p = ', num2str(CC(2,subarray(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    end
%     if ii == size(boxSMP,2)
%         snowLabs = {'Rut','Belly','Virgin Snow'};
%         subplot(5,4,[12,20])
%         boxplot(boxRamm,boxGroup(:,2),'notch','on','Labels',snowLabs)
% %         xtickangle(-10)
%         ylabel({'Strength Index'},'rotation',270,'units','normalized','position',[1.25,.5,0])
%         xlabel('Snow Condition','units','normalized','position',[.5,-.25,0])
%         set(gca,'YAxisLocation','right','fontweight','bold','units','normalized','position',[0.75,0.1125,0.155,0.44])
%         ylim([0,150])
% %         suplabel('SMP','x')
% %         suplabel('Shear Vane','y')
%     end
end

 %% Bootstrapping of ShearVane and SMP Belly Data
% for ii = 1:11
%     for jj = 1:250
% %         tmp1 = datasample(Mbelly(:,ii),5,'Replace',false);
%         tmp1 = datasample(Mbelly(:,ii),3,'Replace',false);
%         tmp2 = datasample(Mbelly(:,ii),3,'Replace',false);
% %         tmp2 = datasample(bellyVane,5,'Replace',false);
%         [tmp3,tmp4] = corrcoef(tmp1,tmp2);
% %         [tmp3,tmp4] = corrcoef(tmp1);    
%         bellyCC(jj,ii) = tmp3(2,1); bellyP(jj,ii) = tmp4(2,1);
%     end
% end
% bellyCC = cell(2,1); bellyP = cell(2,1);
% for kk = 1:2
%     for ii = 1:11
%         for jj = 1:250
% %                     tmp1 = datasample(Mbelly{kk}(:,ii),5,'Replace',false);
%             tmp1 = Mbelly{kk}(:,ii);
% %             tmp2 = datasample(Mbelly{kk}(:,ii),3,'Replace',false);
%                     tmp2 = datasample(shearVaneBellyCelly{kk},5,'Replace',false);
%             [tmp3,tmp4] = corrcoef(tmp1,tmp2);
%             %         [tmp3,tmp4] = corrcoef(tmp1);
%             bellyCC{kk}(jj,ii) = tmp3(2,1); bellyP{kk}(jj,ii) = tmp4(2,1);
%         end
%     end
% end

% 
% for kk = 1:2
%     [~, medVaneIx] = sort(abs(shearVaneBellyCelly{kk} - median(shearVaneBellyCelly{kk})));
%     medVane = shearVaneBellyCelly{kk}(medVaneIx(1:size(medM{kk},1)));
%     for ii = 1:11
%     [~, medBellyIx] = sort(abs(medM{kk}(:,ii) - median(medM{kk}(:,ii))));
%     medBelly = medM{kk}(medBellyIx,ii);
%         for jj = 1:250
%             tmp1 = medM{kk}(:,ii);
% %             tmp1 = medBelly;
% %                         tmp2 = medVane;
%             tmp2 = datasample(medVane,5,'Replace',false);
% %             [tmp3,tmp4] = corrcoef(tmp1,tmp2);
%               [tmp3,tmp4] = corr(tmp1,tmp2(:),'Type','Spearman');
% %                         bellyCC{kk}(ii) = tmp3(2,1); bellyP{kk}(ii) = tmp4(2,1);
% %            bellyCC{kk}(ii) = tmp3; bellyP{kk}(ii) = tmp4;
%             bellyCC{kk}(jj,ii) = tmp3; bellyP{kk}(jj,ii) = tmp4;
%         end
%     end
% end

% figure();subplot(1,2,1);boxplot(bellyCC{1},'Labels',invSMP{1}.vars);
% subplot(1,2,2);boxplot(bellyCC{2},'Labels',invSMP{1}.vars);