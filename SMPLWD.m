%% SMP Comparison with Lightweight Deflectometer
% Tate Meehan 6/5/18
clear; close all; clc;
isWindows = 1;
isLinux = 0;
%% Get LWD Data
isImport = 1;
if isImport
    load('25LWDallsites.mat');
    LWDfiles = {LWD.files}; LWDlabels = {LWD.labels}; stress= {LWD.stress};
        strain= {LWD.strain}; E = cat(2,LWD.E); c = cat(2,LWD.c);
        CC = cat(3,LWD.correlation);
    LWDlocs = 1:length(LWD);
    groups = repmat(LWDlocs,250,1); groups = groups(:);
else
LWDdir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\from_CRREL\NATC 18 West Yellowstone\LWD\lwddata';
LWDfilename = 'ISM LWD NATC 25JAN2018 Compacted Snow Runway Taxiway 05092018.xls';
LWDlocs = 1:10;
LWDnames = [9,25,41,58,74,90,106,122,138,154];
% LWDrows = [{9:24};{25:40};{41:57};{58:73};{74:89};{90:105};{106:121};{122:137};{138:153};{154:169}];
LWDrows = [9,24;25,40;41,57;58,73;74,89;90,105;106,121;122,137;138,153;154,169];

% Remove sites 6,7,8. E bad modulus 
% LWDlocs = 1:7;
% LWDnames = [9,25,41,58,74,138,154];
% % LWDrows = [{9:24};{25:40};{41:57};{58:73};{74:89};{90:105};{106:121};{122:137};{138:153};{154:169}];
% LWDrows = [9,24;25,40;41,57;58,73;74,89;138,153;154,169];
groups = [];
CC = zeros(250,4,length(LWDlocs));
for ii = LWDlocs
[~,LWDfiles(ii)] = xlsread(fullfile(LWDdir,LWDfilename),'NATC Def Corrected',['B',num2str(LWDnames(ii))]);
stress{ii} = xlsread(fullfile(LWDdir,LWDfilename),'NATC Def Corrected',['K',num2str(LWDrows(ii,1)),':','K',num2str(LWDrows(ii,2))]);
stress{ii}(isnan(stress{ii})) = [];
strain{ii} = xlsread(fullfile(LWDdir,LWDfilename),'NATC Def Corrected',['L',num2str(LWDrows(ii,1)),':','L',num2str(LWDrows(ii,2))]);
strain{ii}(isnan(strain{ii})) = [];
% Estimate Young's Modulus
% G = md

for jj = 1:250
%     tmp = -1;
%     maxiter = 1;
%     while tmp<0 && maxiter < 50
        bootstrap = 1:length(strain{ii});
        bootIx = datasample(bootstrap,2,'Replace',false);
        bootstrap(bootIx) = [];
        G = [ones(length(strain{ii}(bootstrap)),1),strain{ii}(bootstrap)];
        d = stress{ii}(bootstrap);
        m = G\d;
        E(jj,ii) = m(2);
        c(jj,ii) = m(1);
        tmp = m(2);
        groups = [groups;ii];
        % Correlation Coefficient
        [R,PP,RL,RU] = corrcoef(strain{ii}(bootstrap),stress{ii}(bootstrap));
        
        CC(jj,1,ii) = R(1,2);
        CC(jj,2,ii) = PP(1,2);
        CC(jj,3,ii) = RL(1,2);
        CC(jj,4,ii) = RU(1,2);
%         maxiter = maxiter+1;
%     end

% Try IRLS Inversion
% G = [ones(length(strain{ii}(bootstrap)),1),strain{ii}(bootstrap)];
% d = stress{ii}(bootstrap);
% x = irls(G,d,1.5E-1,5E-2,1,50);
% Eo(jj,ii) = x(2);
% co(jj,ii) = x(1);

% This is Dope
% mdl = fitlm(strain{ii}(bootstrap),d); % not robust
% mdlr = fitlm(strain{ii}(bootstrap),d,'RobustOpts','on');
end
end
end
%% Visualize
% LWDlabels = [LWDfiles{1},LWDfiles{2},LWDfiles{3},LWDfiles{4},LWDfiles{5},LWDfiles{6},LWDfiles{7}];
% LWDlabels = cell2mat(LWDfiles{:});
% LWDlabels = cat(3,LWDfiles{:});
for ii = 1:length(LWDfiles)
LWDlabels{ii} = LWDfiles{ii}(7:end);
end

figure();
if length(LWDlocs) == 7
for ii = 1:7
    subplot(3,3,ii)
    plot(strain{ii},stress{ii},'ok');hold on;
    medE = quantile(E(:,ii),0.5);
    medC = quantile(c(:,ii),0.5);
%     medEo = quantile(Eo(:,ii),0.5);
%     medCo = quantile(co(:,ii),0.5);
    strainAx = [min(strain{ii}),max(strain{ii})];
    stressL2 = medC + strainAx.*medE;
    plot(strainAx,stressL2,'k');
        title(LWDfiles{ii})
    text(0.75,0.375,['E^{*} = ', num2str(quantile(E(:,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')    
    text(0.75,0.275,['R = ', num2str(quantile(CC(:,1,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.75,0.175,['p = ', num2str(quantile(CC(:,2,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')
%     stressL1=  medCo + strainAx.*medEo;
%     plot(strainAx,stressL1,'--k');
    if ii == 7
        subplot(3,3,[8,9])
        boxplot(E(:),groups,'Labels',LWDlabels)
        xtickangle(-10)
%         subplot(3,3,9)
%         boxplot(Eo(:),groups,'Labels',LWDlabels)
%         xtickangle(-45)
        ylabel('E^{*} (kPa / \mu m)','rotation',270,'units','normalized','position',[1.075,.5,0])
        xlabel('Location','units','normalized','position',[.5,-.25,0])
        set(gca,'fontweight','bold')
    end
    
end
end
if length(LWDlocs) == 10
for ii = 1:10
    subplot(4,3,ii)
    plot(strain{ii},stress{ii},'ok');hold on;
    medE(ii) = quantile(E(:,ii),0.5);
    medC(ii) = quantile(c(:,ii),0.5);
    medP(ii) = quantile(CC(:,2,ii),.5);
%     medEo = quantile(Eo(:,ii),0.5);
%     medCo = quantile(co(:,ii),0.5);
    strainAx = [min(strain{ii}),max(strain{ii})];
    stressL2 = medC(ii) + strainAx.*medE(ii);
    plot(strainAx,stressL2,'k');
%     stressL1=  medCo + strainAx.*medEo;
%     plot(strainAx,stressL1,'--k');
        title(LWDfiles{ii})
    text(0.75,0.475,['E^{*} = ', num2str(quantile(E(:,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')    
    % ,'\plus',num2str(quantile(E(:,ii),.95)),'\minus',num2str(quantile(E(:,ii),.05))
    text(0.75,0.325,['R = ', num2str(quantile(CC(:,1,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.75,0.175,['p = ', num2str(quantile(CC(:,2,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')
    set(gca,'fontweight','bold')
%     axis equal
    if ii == 10
        ylabel('Stress (kPa)');xlabel('Deflection (\mu m)')
        subplot(4,3,[11,12])
        boxplot(E(:),groups,'Labels',LWDlabels); hold on
        plot(0:15,zeros(16,1),'--k')
        xtickangle(-10)
%         subplot(3,3,9)
%         boxplot(Eo(:),groups,'Labels',LWDlabels)
%         xtickangle(-45)
            set(gca,'fontweight','bold')
            ylabel('E^{*} (kPa / \mu m)','rotation',270,'units','normalized','position',[1.075,.5,0])
        xlabel('Location','units','normalized','position',[.5,-.25,0])
        set(gca,'YAxisLocation','right','fontweight','bold')
    end

end
end

%% Export .mat
isExport = 0;
if isExport
LWD = struct('files',LWDfiles,'labels',LWDlabels,'stress',stress,'strain',strain,'E',E,'c',c,'correlation',CC);
for ii = LWDlocs
    LWD(ii).correlation = squeeze(LWD(ii).correlation(:,:,ii));
    LWD(ii).E = squeeze(LWD(ii).E(:,ii));
    LWD(ii).c = squeeze(LWD(ii).c(:,ii));
end
    save('25LWDallsites.mat','LWD','-V7.3')
end
    
%% Load SMP Data and Inversion Result
locs = [1,2,3,4,5];
LWDlocs = find(medP < 0.05 & medE>0); % LWD locations with Significant Correlation at the 5% level
clear medE
 workingDir = pwd;
 medM = cell(length(locs),1);
 nn = 0;
%  shearVane = [];
SVtag = cell(1,1);
q = 0.5; % Test Quantile for Correlations
depth = 50; % Depth for Correlations
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
        SMPname = char({files.name});
        for ii = 1:length(SMP)
            good(ii) = strcmp(SMPname(ii,10:12),'100');
        end
        rutFiles = find(good);
        
        % Shear Vane Depth
%         depth = 50;
        % Runway Ruts
        nn = nn + 1;
        % Rut files
        %         rutVane = [73,93,80,51,52,88,83,98,92];
        %         rutFiles = [45:64];
        badIx = find(rutFiles == 51);%rmv 51
        rutFiles(badIx) = [];
        rutIx = rutFiles;
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
            %             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx),0.34);
            %             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx),0.68);
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
        SMPname = char({files.name});
        for ii = 1:length(SMP)
            good(ii) = strcmp(SMPname(ii,15),'W');
        end
        rutFiles = find(good);
        
        % Shear Vane Depth
%         depth = 50;
        % Runway Ruts
        nn = nn + 1;
        % Rut files
        %         rutVane = [73,93,80,51,52,88,83,98,92];
        %         rutFiles = [45:64];
        badIx = find(rutFiles == 51);%rmv 51
        rutFiles(badIx) = [];
        rutIx = rutFiles;
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
            %             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx),0.34);
            %             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx),0.68);
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
    
    if ff == 3
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
        SMPname = char({files.name});
        for ii = 1:length(SMP)
            good(ii) = strcmp(SMPname(ii,15),'E');
        end
        rutFiles = find(good);
        
        % Shear Vane Depth
%         depth = 50;
        % Runway Ruts
        nn = nn + 1;
        % Rut files
        %         rutVane = [73,93,80,51,52,88,83,98,92];
        %         rutFiles = [45:64];
        badIx = find(rutFiles == 51);%rmv 51
        rutFiles(badIx) = [];
        rutIx = rutFiles;
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
            %             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx),0.34);
            %             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx),0.68);
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
    
    if ff == 4
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
        SMPname = char({files.name});
        for ii = 1:length(SMP)
            good(ii) = strcmp(SMPname(ii,15),'E');
        end
        rutFiles = find(good);
        
        % Shear Vane Depth
%         depth = 50;
        % Runway Ruts
        nn = nn + 1;
        % Rut files
        %         rutVane = [72,48,6,130,55];
        rutFiles = [16:30];
        %         badIx = find(rutFiles == 51);%rmv 51
        %         rutFiles(badIx) = [];
        rutIx = rutFiles;
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
            %             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx),0.34);
            %             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx),0.68);
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
        SMPname = char({files.name});
        for ii = 1:length(SMP)
            good(ii) = strcmp(SMPname(ii,15),'W');
        end
        rutFiles = find(good);
        
        % Shear Vane Depth
%         depth = 50;
        % Runway Ruts
        nn = nn + 1;
        % Rut files
        %         rutVane = [72,48,6,130,55];
        rutFiles = [16:30];
        %         badIx = find(rutFiles == 51);%rmv 51
        %         rutFiles(badIx) = [];
        rutIx = rutFiles;
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
            %             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx),0.34);
            %             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx),0.68);
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
%% Median Statistic of Shear Vane and SMP for Correlation
boxLWD = [];
boxLWDIx = [];
boxSMP = [];
boxSMPIx = [];
for ii = 1:size(medM,1)
    medSMP(ii,:) = nanmedian(medM{ii}(:,:,1),1);
    SMP5(ii,:) = nanmedian(medM{ii}(:,:,2),1);
    SMP95(ii,:) = nanmedian(medM{ii}(:,:,3),1);
    medE(ii,1) = median(E(:,LWDlocs(ii)));
    E5(ii,1) = quantile(E(:,LWDlocs(ii)),0.05);
    E95(ii,1) = quantile(E(:,LWDlocs(ii)),0.95);
    boxLWD = [boxLWD;E(:,LWDlocs(ii))];
    boxLWDIx = [boxLWDIx;LWDlocs(ii)*ones(length(E(:,LWDlocs(ii))),1)];
    boxSMP = [boxSMP;medM{ii}];
    boxSMPIx = [boxSMPIx;ii*ones(size(medM{ii},1),1)];
end
% Compute Correlation Coefficient
for ii = 1:size(medSMP,2)
    [R,PP,RL,RU] = corrcoef(medSMP(:,ii),medE);
    XC(1,ii) = R(1,2);
    XC(2,ii) = PP(1,2);
    XC(3,ii) = RL(1,2);
    XC(4,ii) = RU(1,2);
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
titleVar = cat(2,[invSMP{1}.vars],[invSMP{1}.vars2],'F');
titleVar{8} = 'N_{T}';titleVar{12} = '\delta';titleVar{14} = '\rho';
titleVar{10} = 'F_{m}'; titleVar{15} = 'T_{I}';titleVar{11} = 'F_{med}';
titlevar{9} = 'N_{a}';titleVar{3} = 'P_{c1}';titleVar{4} = 'P_{c2}';
titleVar{13} = 'N_e';titleVar{9} = 'N_a';titleVar{16} = 'N_{m}';
titleStr = {['Characteristic Element Length'],['Rupture Force per mm'],['Probability of Contact'],['Probability of Contact'],['Stiffness'],...
    ['Elastic Modulus'],['Compressive Strength'], ['Number of Ruptures per Characteristic Length'],['Number of Available Elements'],['Mean Rupture Force per mm'],...
    ['Median Rupture Force per mm'],['Deflection'],['Number of Engaged Elements'],['Density'],['Textural Index'],['Number of Measured Ruptures per mm'],['Penetration Force']};
plotarray = [1,2,3,4,5,6,7,8,9,10,11,13,14,15,17,18,19];
subarray = [1,2,12,3,4,13,9,17,10,11,14,15,8,16,5,6,7];
for ii = 1:size(boxSMP,2)
    subplot(5,4,plotarray(ii))
%     plot(medSMP(:,ii),medE,'ok')
    plot(medSMP(:,subarray(ii)),medE,'ok')
    hold on; lsline
    set(gca,'fontweight','bold')
%     title(titleVar{ii})
    title(titleVar{subarray(ii)})

%     text(0.675,0.875,['R = ', num2str(CC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
% %     text(0.675,0.675,['p = ', num2str(CC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.675,0.375,['R = ', num2str(XC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.675,0.175,['p = ', num2str(XC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')

    text(0.675,0.375,['R = ', num2str(XC(1,subarray(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.675,0.175,['p = ', num2str(XC(2,subarray(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    
    if ii == size(boxSMP,2)
        subplot(5,4,[12,20])
        boxplot(boxLWD,boxLWDIx,'Labels',char(LWDlabels{LWDlocs}),'notch','on')
        xtickangle(-20)
        ylabel('E^{*} (kPa / \mu m)','rotation',270,'units','normalized','position',[1.2,.5,0])
        xlabel('Location','units','normalized','position',[.5,-.25,0])
        set(gca,'YAxisLocation','right','fontweight','bold')
%         ylim([0.1,0.45])
%         suplabel('SMP','x')
%         suplabel('Shear Vane','y')
    end
end

% % Significant Correlations
% sigIx = find(CC(2,:) <0.1);
% if length(sigIx)<4;
%     subcols = length(sigIx);
% else
%     subcols  = 4;
% end
% subrows = ceil(length(sigIx)/subcols);
% 
% figure();
% for ii = 1:length(sigIx)
%     subplot(subrows + 1,subcols,ii)
%     plot(medSMP(:,sigIx(ii)),medE,'ok'); hold on;
%     lsline;
% %     plot(SMP5(:,sigIx(ii)),medShearVane,'or','markerfacecolor','r');
% %     plot(SMP95(:,sigIx(ii)),medShearVane,'or','markerfacecolor','r');
%     errorbar(medSMP(:,sigIx(ii)),medE,(medE-E5),-(medE-E95),(medSMP(:,sigIx(ii))-SMP5(:,sigIx(ii))),-(medSMP(:,sigIx(ii))-SMP95(:,sigIx(ii))),'ok')
% %     errorbar(medSMP(:,sigIx(ii)),medShearVane,(medSMP(:,sigIx(ii))-SMP5(:,sigIx(ii))),-(medSMP(:,sigIx(ii))-SMP95(:,sigIx(ii))),'horizontal','ok')
% 
%     text(0.775,0.925,['R = ', num2str(CC(1,sigIx(ii)),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.775,0.85,['p = ', num2str(CC(2,sigIx(ii)),'%4.2f')],'units','normalized','fontweight','bold')
%     title(titleVar{sigIx(ii)})
%     set(gca,'fontweight','bold')
%     ylim([0,135])
%     ylabel('Shear Vane (kPa)')
%     xlabel(titleStr{sigIx(ii)})
%     
% %     axis square
%     % add a box plot
%         if ii == length(sigIx)
%         subplot(subrows + 1,subcols,[ii+1:ii+subcols])
%         boxplot(boxLWD,boxLWDIx,'Labels',SVtag)
%         xtickangle(-10)
%         ylabel('Shear Vane (kPa)','rotation',270,'units','normalized','position',[1.0625,.5,0])
%         xlabel('Location','units','normalized','position',[.5,-.1,0])
%         set(gca,'YAxisLocation','right','fontweight','bold')
% %         suplabel('SMP','x')
% %         suplabel('Shear Vane','y')
%         end
% end
