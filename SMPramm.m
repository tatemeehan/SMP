%% SMP Comparison with Ramm Sonde
% Tate Meehan 6/5/18
clear; close all; clc;
isWindows = 1;
isLinux = 0;
%% Get ramm Data
isImport = 1;
if isImport
    load('rammallsites.mat');
    rammfiles = {ramm.files}; rammnote = {ramm.note}; rammdepth = {ramm.penetration};
        rammindex = {ramm.index}; 
    rammlocs = 1:length(rammfiles);
else
rammdir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\from_CRREL\NATC 18 West Yellowstone\Rammsonde';
rammfilename = 'Ramm NATC MI Feb 2018 TM.xlsx';
sheet = 'ramm';

rammnames = [2,6,16,19,30,33,45,48,56,59,64,73,84,92,102,109,116,126,...
    140,150,156,167,174,188,199,206,225,248,266,272,287,293,308,322,330,...
    336,348,361,368,382,385,390,400,403,410,415,419,428,436,446];
% rammrows = [{9:24};{25:40};{41:57};{58:73};{74:89};{90:105};{106:121};{122:137};{138:153};{154:169}];
rows = [3,14,17,28,31,43,46,54,57,62,71,82,90,100,107,114,124,135,148,154,165,...
    172,186,197,204,223,246,264,270,285,291,306,320,328,334,343,359,366,377,383,388,...
    392,401,408,413,417,426,434,444,452];
rammrows = [rammnames',rows'];
rammlocs = 1:length(rammnames);

% Remove sites 6,7,8. E bad modulus 
% rammlocs = 1:7;
% rammnames = [9,25,41,58,74,138,154];
% % rammrows = [{9:24};{25:40};{41:57};{58:73};{74:89};{90:105};{106:121};{122:137};{138:153};{154:169}];
% rammrows = [9,24;25,40;41,57;58,73;74,89;138,153;154,169];
% groups = [];
% CC = zeros(250,4,length(rammlocs));
for ii = rammlocs
[~,~,tmpfiles] = (xlsread(fullfile(rammdir,rammfilename),sheet,['B',num2str(rammnames(ii)),':','C',num2str(rammnames(ii))]));
%# find numbers
%# convert to string
tmpfiles{1} = [tmpfiles{1},num2str(tmpfiles{2})];
rammfiles{ii} = tmpfiles{1};
rammindex{ii} = xlsread(fullfile(rammdir,rammfilename),sheet,['J',num2str(rammrows(ii,1)+1),':','J',num2str(rammrows(ii,2))]);
rammstrength{ii} = 4.78*log(rammindex{ii})-14.72;
indexix = find(rammindex{ii});
rammindex{ii}=rammindex{ii}(indexix);
rammdepth{ii} = xlsread(fullfile(rammdir,rammfilename),sheet,['F',num2str(rammrows(ii,1)+1),':','F',num2str(rammrows(ii,2))]);
rammdepth{ii} = rammdepth{ii}(indexix);
[~,rammnote{ii}] = xlsread(fullfile(rammdir,rammfilename),sheet,['C',num2str(rammnames(ii)+1)]);
medrammindex{ii} = quantile(rammindex{ii},.5);
maxrammindex{ii} = quantile(rammindex{ii},.995);
minrammindex{ii} = quantile(rammindex{ii},.005);
clear indexix

% Estimate Young's Modulus
% G = md

% for jj = 1:250
% %     tmp = -1;
% %     maxiter = 1;
% %     while tmp<0 && maxiter < 50
%         bootstrap = 1:length(strain{ii});
%         bootIx = datasample(bootstrap,2,'Replace',false);
%         bootstrap(bootIx) = [];
%         G = [ones(length(strain{ii}(bootstrap)),1),strain{ii}(bootstrap)];
%         d = stress{ii}(bootstrap);
%         m = G\d;
%         E(jj,ii) = m(2);
%         c(jj,ii) = m(1);
%         tmp = m(2);
%         groups = [groups;ii];
%         % Correlation Coefficient
%         [R,PP,RL,RU] = corrcoef(strain{ii}(bootstrap),stress{ii}(bootstrap));
%         
%         CC(jj,1,ii) = R(1,2);
%         CC(jj,2,ii) = PP(1,2);
%         CC(jj,3,ii) = RL(1,2);
%         CC(jj,4,ii) = RU(1,2);
% %         maxiter = maxiter+1;
% %     end
% 
% % Try IRLS Inversion
% % G = [ones(length(strain{ii}(bootstrap)),1),strain{ii}(bootstrap)];
% % d = stress{ii}(bootstrap);
% % x = irls(G,d,1.5E-1,5E-2,1,50);
% % Eo(jj,ii) = x(2);
% % co(jj,ii) = x(1);
% 
% % This is Dope
% % mdl = fitlm(strain{ii}(bootstrap),d); % not robust
% % mdlr = fitlm(strain{ii}(bootstrap),d,'RobustOpts','on');
% end
end
end
%% Visualize
% % rammlabels = [rammfiles{1},rammfiles{2},rammfiles{3},rammfiles{4},rammfiles{5},rammfiles{6},rammfiles{7}];
% % rammlabels = cell2mat(rammfiles{:});
% % rammlabels = cat(3,rammfiles{:});
% for ii = 1:length(tmpfiles)
% rammlabels{ii} = tmpfiles{ii}(7:end);
% end
% 
% figure();
% if length(rammlocs) == 7
% for ii = 1:7
%     subplot(3,3,ii)
%     plot(strain{ii},stress{ii},'ok');hold on;
%     medE = quantile(E(:,ii),0.5);
%     medC = quantile(c(:,ii),0.5);
% %     medEo = quantile(Eo(:,ii),0.5);
% %     medCo = quantile(co(:,ii),0.5);
%     strainAx = [min(strain{ii}),max(strain{ii})];
%     stressL2 = medC + strainAx.*medE;
%     plot(strainAx,stressL2,'k');
%         title(tmpfiles{ii})
%     text(0.75,0.375,['E^{*} = ', num2str(quantile(E(:,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')    
%     text(0.75,0.275,['R = ', num2str(quantile(CC(:,1,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.75,0.175,['p = ', num2str(quantile(CC(:,2,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')
% %     stressL1=  medCo + strainAx.*medEo;
% %     plot(strainAx,stressL1,'--k');
%     if ii == 7
%         subplot(3,3,[8,9])
%         boxplot(E(:),groups,'Labels',rammlabels)
%         xtickangle(-10)
% %         subplot(3,3,9)
% %         boxplot(Eo(:),groups,'Labels',rammlabels)
% %         xtickangle(-45)
%         ylabel('E^{*} (kPa / \mu m)','rotation',270,'units','normalized','position',[1.075,.5,0])
%         xlabel('Location','units','normalized','position',[.5,-.25,0])
%         set(gca,'fontweight','bold')
%     end
%     
% end
% end
% if length(rammlocs) == 10
% for ii = 1:10
%     subplot(4,3,ii)
%     plot(strain{ii},stress{ii},'ok');hold on;
%     medE(ii) = quantile(E(:,ii),0.5);
%     medC(ii) = quantile(c(:,ii),0.5);
%     medP(ii) = quantile(CC(:,2,ii),.5);
% %     medEo = quantile(Eo(:,ii),0.5);
% %     medCo = quantile(co(:,ii),0.5);
%     strainAx = [min(strain{ii}),max(strain{ii})];
%     stressL2 = medC(ii) + strainAx.*medE(ii);
%     plot(strainAx,stressL2,'k');
% %     stressL1=  medCo + strainAx.*medEo;
% %     plot(strainAx,stressL1,'--k');
%         title(tmpfiles{ii})
%     text(0.75,0.475,['E^{*} = ', num2str(quantile(E(:,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')    
%     % ,'\plus',num2str(quantile(E(:,ii),.95)),'\minus',num2str(quantile(E(:,ii),.05))
%     text(0.75,0.325,['R = ', num2str(quantile(CC(:,1,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.75,0.175,['p = ', num2str(quantile(CC(:,2,ii),.5),'%4.2f')],'units','normalized','fontweight','bold')
%     set(gca,'fontweight','bold')
% %     axis equal
%     if ii == 10
%         ylabel('Stress (kPa)');xlabel('Deflection (\mu m)')
%         subplot(4,3,[11,12])
%         boxplot(E(:),groups,'Labels',rammlabels); hold on
%         plot(0:15,zeros(16,1),'--k')
%         xtickangle(-10)
% %         subplot(3,3,9)
% %         boxplot(Eo(:),groups,'Labels',rammlabels)
% %         xtickangle(-45)
%             set(gca,'fontweight','bold')
%             ylabel('E^{*} (kPa / \mu m)','rotation',270,'units','normalized','position',[1.075,.5,0])
%         xlabel('Location','units','normalized','position',[.5,-.25,0])
%         set(gca,'YAxisLocation','right','fontweight','bold')
%     end
% 
% end
% end

%% Export .mat
isExport = 0;
if isExport
ramm = struct('files',rammfiles,'note',rammnote,'penetration',rammdepth,'index',rammindex,'medindex',medrammindex,'minindex',minrammindex,'maxindex',maxrammindex, 'strength',rammstrength);

    save('rammallsites.mat','ramm','-V7.3')
end
    
%% Load SMP Data and Inversion Result
locs = [1,2,3,4,5,6];
boxRamm = [];
boxGroup = [];
boxRstrength = [];
 workingDir = pwd;
 medM = cell(length(locs),1);
 normMedM = cell(length(locs),1);
 subNormMedM = cell(length(locs),1);
 nn = 0;
%  shearVane = [];
SVtag = cell(1,1);
q = 0.5; % Test Quantile for Correlations
depth = 50; % Depth for Correlations
for ff = locs
    % 27Jan2018 data
    if ff == 1
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/27Jan2018/SMPS/25CRWYSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\27Jan2018\SMPS\27SWVSSPMS';
        end
        
        tag = '27SWVSSPMS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        SVtag{1} = [SVtag{1},{'27SWVS'}];
        SMPname = char({files.name});
%         for ii = 1:length(SMP)
%             good(ii) = strcmp(SMPname(ii,10:12),'100');
%         end
%         rutFiles = find(good);
        % Ramm Files
        rammrut = [2,4];
        rammbelly = [3];
        rammvs = [1,5];
        
        % SMP Files
        smprut = [274,275,276,277,278,284,285,286,287,288,289];
        smpbelly = [279,280,281,282,283];
        smpvs = [269,270,271,272,273,290,291,292,293,294];
        smpnums = str2num(SMPname(:,6:8));
        for ii = 1:length(smprut)
            try
            smprutix(ii) = find(smprut(ii) == smpnums);
            catch
            end
        end
        smprutix(smprutix==0) = [];
        for ii = 1:length(smpbelly)
            try
            smpbellyix(ii) = find(smpbelly(ii) == smpnums);
            catch
            end
        end   
        smpbellyix(smpbellyix==0) = [];
        
        for ii = 1:length(smpvs)
            try
            smpvsix(ii) = find(smpvs(ii) == smpnums);
            catch
            end
        end 
        smpvsix(smpvsix==0) = []; 
        
        % Ramm Sonde Depth
        for ii = 1:length(rammrut)
            depthrut(ii) = median(ramm(rammrut(ii)).penetration);
            medrammrut(ii) = median(ramm(rammrut(ii)).index);
            boxRamm = [boxRamm;ramm(rammrut(ii)).index];
            medstrengthrut(ii) = median(ramm(rammrut(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammrut(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammrut(ii)).index),1)...
                ,1*ones(length(ramm(rammrut(ii)).index),1)]];
        end
        depthrut = median(depthrut);
        
        for ii = 1:length(rammbelly)
            depthbelly(ii) = median(ramm(rammbelly(ii)).penetration);
            medrammbelly(ii) = median(ramm(rammbelly(ii)).index);
            boxRamm = [boxRamm;ramm(rammbelly(ii)).index];
            medstrengthbelly(ii) = median(ramm(rammbelly(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammbelly(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammbelly(ii)).index),1)...
                ,2*ones(length(ramm(rammbelly(ii)).index),1)]];            
        end
        depthbelly = median(depthbelly);
        for ii = 1:length(rammvs)
            depthvs(ii) = median(ramm(rammvs(ii)).penetration);
            medrammvs(ii) = median(ramm(rammvs(ii)).index);
            boxRamm = [boxRamm;ramm(rammvs(ii)).index];
            medstrengthvs(ii) = median(ramm(rammvs(ii)).strength);            
            boxRstrength = [boxRstrength;ramm(rammvs(ii)).strength];            
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammvs(ii)).index),1)...
                ,3*ones(length(ramm(rammvs(ii)).index),1)]];            
        end
        depthvs = median(depthvs);

        nn = nn + 1;
        % Loop over Rut Files
        smpIx = smprutix;
        depth = depthrut*10;
        quantileSMP
        medramm(nn) = median(medrammrut);
        medstrength(nn) = median(medstrengthrut);
        
        nn = nn + 1;
        % Loop Over Belly Files
        smpIx = smpbellyix;
        depth = depthbelly*10;
        quantileSMP
        medramm(nn) = median(medrammbelly);
        medstrength(nn) = median(medstrengthbelly);
        
        nn = nn + 1;
        % Loop Over VS Files
        smpIx = smpvsix;
        depth = depthvs*10;
        quantileSMP
        medramm(nn) = median(medrammvs);
        medstrength(nn) = median(medstrengthvs);
        
        clear ('smprutix','smpbellyix','smpvsix');
        
%         nn = nn + 1;        
    end
    
    if ff == 2
        % 28Jan2018 data
        % 28CMobLPSMPS
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/25Jan2018/SMPS/25CRWYSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\28Jan2018\SMPS\28CMobLPSMPS';
        end
        
        tag = '28CMobLPSMPS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        SVtag{1} = [SVtag{1},{'28VSLVSR'}];
        SMPname = char({files.name});

        
        % 28VSLVSR
        % Ramm Files
        rammrut = [8,11];
        rammbelly = [7,9,10];
        rammvs = [6];
        
        % SMP Files
        smprut = [337,338,339,340,341,342,343,344,345,351,352,353,354,355];
        smpbelly = [332,333,334,335,336,346,347,348,349,350,356,357,358,359,360];
        smpvs = [327,328,329,330,331];
        smpnums = str2num(SMPname(:,6:8));
        for ii = 1:length(smprut)
            try
            smprutix(ii) = find(smprut(ii) == smpnums);
            catch
            end
        end
        smprutix(smprutix==0) = [];
        for ii = 1:length(smpbelly)
            try
            smpbellyix(ii) = find(smpbelly(ii) == smpnums);
            catch
            end
        end   
        smpbellyix(smpbellyix==0) = [];
        
        for ii = 1:length(smpvs)
            try
            smpvsix(ii) = find(smpvs(ii) == smpnums);
            catch
            end
        end 
        smpvsix(smpvsix==0) = [];
        
        % Ramm Sonde Depth
        for ii = 1:length(rammrut)
            depthrut(ii) = median(ramm(rammrut(ii)).penetration);
            medrammrut(ii) = median(ramm(rammrut(ii)).index);
            boxRamm = [boxRamm;ramm(rammrut(ii)).index];
            medstrengthrut(ii) = median(ramm(rammrut(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammrut(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammrut(ii)).index),1)...
                ,1*ones(length(ramm(rammrut(ii)).index),1)]];
        end
        depthrut = median(depthrut);
        
        for ii = 1:length(rammbelly)
            depthbelly(ii) = median(ramm(rammbelly(ii)).penetration);
            medrammbelly(ii) = median(ramm(rammbelly(ii)).index);
            boxRamm = [boxRamm;ramm(rammbelly(ii)).index];
            medstrengthbelly(ii) = median(ramm(rammbelly(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammbelly(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammbelly(ii)).index),1)...
                ,2*ones(length(ramm(rammbelly(ii)).index),1)]];            
        end
        depthbelly = median(depthbelly);
        for ii = 1:length(rammvs)
            depthvs(ii) = median(ramm(rammvs(ii)).penetration);
            medrammvs(ii) = median(ramm(rammvs(ii)).index);
            boxRamm = [boxRamm;ramm(rammvs(ii)).index];
            medstrengthvs(ii) = median(ramm(rammvs(ii)).strength);            
            boxRstrength = [boxRstrength;ramm(rammvs(ii)).strength];            
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammvs(ii)).index),1)...
                ,3*ones(length(ramm(rammvs(ii)).index),1)]];            
        end
        depthvs = median(depthvs);
        
        nn = nn + 1;
        % Loop over Rut Files
        smpIx = smprutix;
        depth = depthrut*10;
        quantileSMP
        medramm(nn) = median(medrammrut);
        medstrength(nn) = median(medstrengthrut);
        
        nn = nn + 1;
        % Loop Over Belly Files
        smpIx = smpbellyix;
        depth = depthbelly*10;
        quantileSMP
        medramm(nn) = median(medrammbelly);
        medstrength(nn) = median(medstrengthbelly);
        
        nn = nn + 1;
        % Loop Over VS Files
        smpIx = smpvsix;
        depth = depthvs*10;
        quantileSMP
        medramm(nn) = median(medrammvs);
        medstrength(nn) = median(medstrengthvs);

        
        clear ('smprutix','smpbellyix','smpvsix');
        
        % 28MobLPLVSR
        SVtag{1} = [SVtag{1},{'28MBLPLVSR'}];
        % Ramm Files
        rammrut = [13,15];
        rammbelly = [12,14];
        rammvs = [16];
        
        % SMP Files
        smprut = [380,381,382,388,389];
        smpbelly = [383,384,385,386,387];
        smpvs = [365,366,367];
        smpnums = str2num(SMPname(:,6:8));
        for ii = 1:length(smprut)
            try
            smprutix(ii) = find(smprut(ii) == smpnums);
            catch
            end
        end
        smprutix(smprutix==0) = [];
        for ii = 1:length(smpbelly)
            try
            smpbellyix(ii) = find(smpbelly(ii) == smpnums);
            catch
            end
        end   
        smpbellyix(smpbellyix==0) = [];
        
        for ii = 1:length(smpvs)
            try
            smpvsix(ii) = find(smpvs(ii) == smpnums);
            catch
            end
        end 
        smpvsix(smpvsix==0) = [];
        
        % Ramm Sonde Depth
        for ii = 1:length(rammrut)
            depthrut(ii) = median(ramm(rammrut(ii)).penetration);
            medrammrut(ii) = median(ramm(rammrut(ii)).index);
            boxRamm = [boxRamm;ramm(rammrut(ii)).index];
            medstrengthrut(ii) = median(ramm(rammrut(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammrut(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammrut(ii)).index),1)...
                ,1*ones(length(ramm(rammrut(ii)).index),1)]];
        end
        depthrut = median(depthrut);
        
        for ii = 1:length(rammbelly)
            depthbelly(ii) = median(ramm(rammbelly(ii)).penetration);
            medrammbelly(ii) = median(ramm(rammbelly(ii)).index);
            boxRamm = [boxRamm;ramm(rammbelly(ii)).index];
            medstrengthbelly(ii) = median(ramm(rammbelly(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammbelly(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammbelly(ii)).index),1)...
                ,2*ones(length(ramm(rammbelly(ii)).index),1)]];            
        end
        depthbelly = median(depthbelly);
        for ii = 1:length(rammvs)
            depthvs(ii) = median(ramm(rammvs(ii)).penetration);
            medrammvs(ii) = median(ramm(rammvs(ii)).index);
            boxRamm = [boxRamm;ramm(rammvs(ii)).index];
            medstrengthvs(ii) = median(ramm(rammvs(ii)).strength);            
            boxRstrength = [boxRstrength;ramm(rammvs(ii)).strength];            
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammvs(ii)).index),1)...
                ,3*ones(length(ramm(rammvs(ii)).index),1)]];            
        end
        depthvs = median(depthvs);
        
        nn = nn + 1;
        % Loop over Rut Files
        smpIx = smprutix;
        depth = depthrut*10;
        quantileSMP
        medramm(nn) = median(medrammrut);
        medstrength(nn) = median(medstrengthrut);
        
        nn = nn + 1;
        % Loop Over Belly Files
        smpIx = smpbellyix;
        depth = depthbelly*10;
        quantileSMP
        medramm(nn) = median(medrammbelly);
        medstrength(nn) = median(medstrengthbelly);
        
        nn = nn + 1;
        % Loop Over VS Files
        smpIx = smpvsix;
        depth = depthvs*10;
        quantileSMP
        medramm(nn) = median(medrammvs);
        medstrength(nn) = median(medstrengthvs);

        
        clear ('smprutix','smpbellyix','smpvsix');
        
        % 28MobLPMRZR
        SVtag{1} = [SVtag{1},{'28MBLPRZR'}];
        % Ramm Files
        rammrut = [18,19];
        rammbelly = [16];
        rammvs = [17];
        
        % SMP Files
        smprut = [369,370,379,390,391,392,398,399];
        smpbelly = [371,373,374,375,376,377,378,393,394,395,396,397];
        smpvs = [365,366,367,368];
        smpnums = str2num(SMPname(:,6:8));
        for ii = 1:length(smprut)
            try
            smprutix(ii) = find(smprut(ii) == smpnums);
            catch
            end
        end
        smprutix(smprutix==0) = [];
        for ii = 1:length(smpbelly)
            try
            smpbellyix(ii) = find(smpbelly(ii) == smpnums);
            catch
            end
        end   
        smpbellyix(smpbellyix==0) = [];
        
        for ii = 1:length(smpvs)
            try
            smpvsix(ii) = find(smpvs(ii) == smpnums);
            catch
            end
        end 
        smpvsix(smpvsix==0) = [];
        % Ramm Sonde Depth
        for ii = 1:length(rammrut)
            depthrut(ii) = median(ramm(rammrut(ii)).penetration);
            medrammrut(ii) = median(ramm(rammrut(ii)).index);
            boxRamm = [boxRamm;ramm(rammrut(ii)).index];
            medstrengthrut(ii) = median(ramm(rammrut(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammrut(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammrut(ii)).index),1)...
                ,1*ones(length(ramm(rammrut(ii)).index),1)]];
        end
        depthrut = median(depthrut);
        
        for ii = 1:length(rammbelly)
            depthbelly(ii) = median(ramm(rammbelly(ii)).penetration);
            medrammbelly(ii) = median(ramm(rammbelly(ii)).index);
            boxRamm = [boxRamm;ramm(rammbelly(ii)).index];
            medstrengthbelly(ii) = median(ramm(rammbelly(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammbelly(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammbelly(ii)).index),1)...
                ,2*ones(length(ramm(rammbelly(ii)).index),1)]];            
        end
        depthbelly = median(depthbelly);
        for ii = 1:length(rammvs)
            depthvs(ii) = median(ramm(rammvs(ii)).penetration);
            medrammvs(ii) = median(ramm(rammvs(ii)).index);
            boxRamm = [boxRamm;ramm(rammvs(ii)).index];
            medstrengthvs(ii) = median(ramm(rammvs(ii)).strength);            
            boxRstrength = [boxRstrength;ramm(rammvs(ii)).strength];            
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammvs(ii)).index),1)...
                ,3*ones(length(ramm(rammvs(ii)).index),1)]];            
        end
        depthvs = median(depthvs);
        
        nn = nn + 1;
        % Loop over Rut Files
        smpIx = smprutix;
        depth = depthrut*10;
        quantileSMP
        medramm(nn) = median(medrammrut);
        medstrength(nn) = median(medstrengthrut);
        
        nn = nn + 1;
        % Loop Over Belly Files
        smpIx = smpbellyix;
        depth = depthbelly*10;
        quantileSMP
        medramm(nn) = median(medrammbelly);
        medstrength(nn) = median(medstrengthbelly);
        
        nn = nn + 1;
        % Loop Over VS Files
        smpIx = smpvsix;
        depth = depthvs*10;
        quantileSMP
        medramm(nn) = median(medrammvs);
        medstrength(nn) = median(medstrengthvs);    
        
        clear ('smprutix','smpbellyix','smpvsix');

    end
    if ff == 3
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
        SMPname = char({files.name});

        % LVSR CD Ruts
        SVtag{1} = [SVtag{1},{'29STXYLVSRCD'}];
        % Ramm Files
        rammrut = [21,25];
        rammbelly = [22,23,24];
        rammvs = [20];
        
        % SMP Files
%                 smprut = [440,449,450];
        smprut = [440,448,449,450];
        smpbelly = [430:434];
        smpvs = [423:426];
        smpnums = str2num(SMPname(:,6:8));
        for ii = 1:length(smprut)
            try
            smprutix(ii) = find(smprut(ii) == smpnums);
            catch
            end
        end
        smprutix(smprutix==0) = [];
        for ii = 1:length(smpbelly)
            try
            smpbellyix(ii) = find(smpbelly(ii) == smpnums);
            catch
            end
        end   
        smpbellyix(smpbellyix==0) = [];
        
        for ii = 1:length(smpvs)
            try
            smpvsix(ii) = find(smpvs(ii) == smpnums);
            catch
            end
        end 
        smpvsix(smpvsix==0) = [];
        % Ramm Sonde Depth
        for ii = 1:length(rammrut)
            depthrut(ii) = median(ramm(rammrut(ii)).penetration);
            medrammrut(ii) = median(ramm(rammrut(ii)).index);
            boxRamm = [boxRamm;ramm(rammrut(ii)).index];
            medstrengthrut(ii) = median(ramm(rammrut(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammrut(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammrut(ii)).index),1)...
                ,1*ones(length(ramm(rammrut(ii)).index),1)]];
        end
        depthrut = median(depthrut);
        
        for ii = 1:length(rammbelly)
            depthbelly(ii) = median(ramm(rammbelly(ii)).penetration);
            medrammbelly(ii) = median(ramm(rammbelly(ii)).index);
            boxRamm = [boxRamm;ramm(rammbelly(ii)).index];
            medstrengthbelly(ii) = median(ramm(rammbelly(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammbelly(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammbelly(ii)).index),1)...
                ,2*ones(length(ramm(rammbelly(ii)).index),1)]];            
        end
        depthbelly = median(depthbelly);
        for ii = 1:length(rammvs)
            depthvs(ii) = median(ramm(rammvs(ii)).penetration);
            medrammvs(ii) = median(ramm(rammvs(ii)).index);
            boxRamm = [boxRamm;ramm(rammvs(ii)).index];
            medstrengthvs(ii) = median(ramm(rammvs(ii)).strength);            
            boxRstrength = [boxRstrength;ramm(rammvs(ii)).strength];            
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammvs(ii)).index),1)...
                ,3*ones(length(ramm(rammvs(ii)).index),1)]];            
        end
        depthvs = median(depthvs);
        
        nn = nn + 1;
        % Loop over Rut Files
        smpIx = smprutix;
        depth = depthrut*10;
        quantileSMP
        medramm(nn) = median(medrammrut);
        medstrength(nn) = median(medstrengthrut);
        
        nn = nn + 1;
        % Loop Over Belly Files
        smpIx = smpbellyix;
        depth = depthbelly*10;
        quantileSMP
        medramm(nn) = median(medrammbelly);
        medstrength(nn) = median(medstrengthbelly);
        
        nn = nn + 1;
        % Loop Over VS Files
        smpIx = smpvsix;
        depth = depthvs*10;
        quantileSMP
        medramm(nn) = median(medrammvs);
        medstrength(nn) = median(medstrengthvs);    
        
        clear ('smprutix','smpbellyix','smpvsix');
    end
    
        % 30Jan2018 data
    if ff == 4
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
        SMPname = char({files.name});

        % Mobility Loop Entrance Ruts
        SVtag{1} = [SVtag{1},{'30NMOBLP'}];
        % Ramm Files
        rammrut = [29,31,34,36];
        rammbelly = [27,28,30,32,35];
        rammvs = [26];
        
        % SMP Files
        smprut = [460,461,467,468,469,486];
        smpbelly = [462,463,464,465,466];
        smpvs = [456,457,458,459];
        smpnums = str2num(SMPname(:,6:8));
        for ii = 1:length(smprut)
            try
            smprutix(ii) = find(smprut(ii) == smpnums);
            catch
            end
        end
        smprutix(smprutix==0) = [];
        for ii = 1:length(smpbelly)
            try
            smpbellyix(ii) = find(smpbelly(ii) == smpnums);
            catch
            end
        end   
        smpbellyix(smpbellyix==0) = [];
        
        for ii = 1:length(smpvs)
            try
            smpvsix(ii) = find(smpvs(ii) == smpnums);
            catch
            end
        end 
        smpvsix(smpvsix==0) = [];
        % Ramm Sonde Depth
        for ii = 1:length(rammrut)
            depthrut(ii) = median(ramm(rammrut(ii)).penetration);
            medrammrut(ii) = median(ramm(rammrut(ii)).index);
            boxRamm = [boxRamm;ramm(rammrut(ii)).index];
            medstrengthrut(ii) = median(ramm(rammrut(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammrut(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammrut(ii)).index),1)...
                ,1*ones(length(ramm(rammrut(ii)).index),1)]];
        end
        depthrut = median(depthrut);
        
        for ii = 1:length(rammbelly)
            depthbelly(ii) = median(ramm(rammbelly(ii)).penetration);
            medrammbelly(ii) = median(ramm(rammbelly(ii)).index);
            boxRamm = [boxRamm;ramm(rammbelly(ii)).index];
            medstrengthbelly(ii) = median(ramm(rammbelly(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammbelly(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammbelly(ii)).index),1)...
                ,2*ones(length(ramm(rammbelly(ii)).index),1)]];            
        end
        depthbelly = median(depthbelly);
        for ii = 1:length(rammvs)
            depthvs(ii) = median(ramm(rammvs(ii)).penetration);
            medrammvs(ii) = median(ramm(rammvs(ii)).index);
            boxRamm = [boxRamm;ramm(rammvs(ii)).index];
            medstrengthvs(ii) = median(ramm(rammvs(ii)).strength);            
            boxRstrength = [boxRstrength;ramm(rammvs(ii)).strength];            
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammvs(ii)).index),1)...
                ,3*ones(length(ramm(rammvs(ii)).index),1)]];            
        end
        depthvs = median(depthvs);
        
        nn = nn + 1;
        % Loop over Rut Files
        smpIx = smprutix;
        depth = depthrut*10;
        quantileSMP
        medramm(nn) = median(medrammrut);
        medstrength(nn) = median(medstrengthrut);
        
        nn = nn + 1;
        % Loop Over Belly Files
        smpIx = smpbellyix;
        depth = depthbelly*10;
        quantileSMP
        medramm(nn) = median(medrammbelly);
        medstrength(nn) = median(medstrengthbelly);
        
        nn = nn + 1;
        % Loop Over VS Files
        smpIx = smpvsix;
        depth = depthvs*10;
        quantileSMP
        medramm(nn) = median(medrammvs);
        medstrength(nn) = median(medstrengthvs);
        
        clear ('smprutix','smpbellyix','smpvsix');
    end
    
        % 30Jan2018 data
    if ff == 5
        % 30NMobLPSPMS data
        tic
        if isLinux
            dataDir = '/home/tatemeehan/CRREL_SnowCompaction/W_YELLOWSTONE/DATA/30Jan2018/SMPS/30STXYMTVRCDSMPS';
        end
        if isWindows
            dataDir = 'E:\CRREL_SnowCompaction\W_YELLOWSTONE\DATA\30Jan2018\SMPS\30STXYMTVRCDSMPS';
        end         
        tag = '30STXYMTVRCDSMPS';
        load(fullfile(dataDir,[tag,'data']));
        load(fullfile(dataDir,[tag,'results']));
        load(fullfile(dataDir,[tag,'files']))
        SVtag{1} = [SVtag{1},{'30SRWYCD'}];
        SMPname = char({files.name});

        % LVSR CD Ruts
        % Ramm Files
        rammrut = [37,39];
        rammbelly = [38];
        rammvs = [40,41];
        
        % SMP Files
        smprut = [508:510];
        smpbelly = [503:507];
        smpvs = [498:502];
        smpnums = str2num(SMPname(:,6:8));
        for ii = 1:length(smprut)
            try
            smprutix(ii) = find(smprut(ii) == smpnums);
            catch
            end
        end
        smprutix(smprutix==0) = [];
        for ii = 1:length(smpbelly)
            try
            smpbellyix(ii) = find(smpbelly(ii) == smpnums);
            catch
            end
        end   
        smpbellyix(smpbellyix==0) = [];
        
        for ii = 1:length(smpvs)
            try
            smpvsix(ii) = find(smpvs(ii) == smpnums);
            catch
            end
        end 
        smpvsix(smpvsix==0) = [];
        % Ramm Sonde Depth
        for ii = 1:length(rammrut)
            depthrut(ii) = median(ramm(rammrut(ii)).penetration);
            medrammrut(ii) = median(ramm(rammrut(ii)).index);
            boxRamm = [boxRamm;ramm(rammrut(ii)).index];
            medstrengthrut(ii) = median(ramm(rammrut(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammrut(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammrut(ii)).index),1)...
                ,1*ones(length(ramm(rammrut(ii)).index),1)]];
        end
        depthrut = median(depthrut);
        
        for ii = 1:length(rammbelly)
            depthbelly(ii) = median(ramm(rammbelly(ii)).penetration);
            medrammbelly(ii) = median(ramm(rammbelly(ii)).index);
            boxRamm = [boxRamm;ramm(rammbelly(ii)).index];
            medstrengthbelly(ii) = median(ramm(rammbelly(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammbelly(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammbelly(ii)).index),1)...
                ,2*ones(length(ramm(rammbelly(ii)).index),1)]];            
        end
        depthbelly = median(depthbelly);
        for ii = 1:length(rammvs)
            depthvs(ii) = median(ramm(rammvs(ii)).penetration);
            medrammvs(ii) = median(ramm(rammvs(ii)).index);
            boxRamm = [boxRamm;ramm(rammvs(ii)).index];
            medstrengthvs(ii) = median(ramm(rammvs(ii)).strength);            
            boxRstrength = [boxRstrength;ramm(rammvs(ii)).strength];            
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammvs(ii)).index),1)...
                ,3*ones(length(ramm(rammvs(ii)).index),1)]];            
        end
        depthvs = median(depthvs);
        
        nn = nn + 1;
        % Loop over Rut Files
        smpIx = smprutix;
        depth = depthrut*10;
        quantileSMP
        medramm(nn) = median(medrammrut);
        medstrength(nn) = median(medstrengthrut);
        
        nn = nn + 1;
        % Loop Over Belly Files
        smpIx = smpbellyix;
        depth = depthbelly*10;
        quantileSMP
        medramm(nn) = median(medrammbelly);
        medstrength(nn) = median(medstrengthbelly);
        
        nn = nn + 1;
        % Loop Over VS Files
        smpIx = smpvsix;
        depth = depthvs*10;
        quantileSMP
        medramm(nn) = median(medrammvs);
        medstrength(nn) = median(medstrengthvs);
        
        clear ('smprutix','smpbellyix','smpvsix');
    end
    
        % 31Jan2018 data
    if ff == 6
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
        SVtag{1} = [SVtag{1},{'31MTVRCD'}];
        SMPname = char({files.name});

        % LVSR CD Ruts
        % Ramm Files
        rammrut = [42,45,46,47,50];
        rammbelly = [43,44,48,49];
        rammvs = [40,41];
        
        % SMP Files
        smprut = [526];%,528];
        smpbelly = [521:525];
        smpvs = [516:520];
        smpnums = str2num(SMPname(:,6:8));
        for ii = 1:length(smprut)
            try
            smprutix(ii) = find(smprut(ii) == smpnums);
            catch
            end
        end
        smprutix(smprutix==0) = [];
        for ii = 1:length(smpbelly)
            try
            smpbellyix(ii) = find(smpbelly(ii) == smpnums);
            catch
            end
        end   
        smpbellyix(smpbellyix==0) = [];
        
        for ii = 1:length(smpvs)
            try
            smpvsix(ii) = find(smpvs(ii) == smpnums);
            catch
            end
        end 
        smpvsix(smpvsix==0) = [];
        % Ramm Sonde Depth
        for ii = 1:length(rammrut)
            depthrut(ii) = median(ramm(rammrut(ii)).penetration);
            medrammrut(ii) = median(ramm(rammrut(ii)).index);
            boxRamm = [boxRamm;ramm(rammrut(ii)).index];
            medstrengthrut(ii) = median(ramm(rammrut(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammrut(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammrut(ii)).index),1)...
                ,1*ones(length(ramm(rammrut(ii)).index),1)]];
        end
        depthrut = median(depthrut);
        
        for ii = 1:length(rammbelly)
            depthbelly(ii) = median(ramm(rammbelly(ii)).penetration);
            medrammbelly(ii) = median(ramm(rammbelly(ii)).index);
            boxRamm = [boxRamm;ramm(rammbelly(ii)).index];
            medstrengthbelly(ii) = median(ramm(rammbelly(ii)).strength);
            boxRstrength = [boxRstrength;ramm(rammbelly(ii)).strength];
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammbelly(ii)).index),1)...
                ,2*ones(length(ramm(rammbelly(ii)).index),1)]];            
        end
        depthbelly = median(depthbelly);
        for ii = 1:length(rammvs)
            depthvs(ii) = median(ramm(rammvs(ii)).penetration);
            medrammvs(ii) = median(ramm(rammvs(ii)).index);
            boxRamm = [boxRamm;ramm(rammvs(ii)).index];
            medstrengthvs(ii) = median(ramm(rammvs(ii)).strength);            
            boxRstrength = [boxRstrength;ramm(rammvs(ii)).strength];            
            boxGroup = [boxGroup;[locs(ff).*ones(length(ramm(rammvs(ii)).index),1)...
                ,3*ones(length(ramm(rammvs(ii)).index),1)]];            
        end
        depthvs = median(depthvs);
        
        nn = nn + 1;
        % Loop over Rut Files
        smpIx = smprutix;
        depth = depthrut*10;
        quantileSMP
        medramm(nn) = median(medrammrut);
        medstrength(nn) = median(medstrengthrut);
        
        nn = nn + 1;
        % Loop Over Belly Files
        smpIx = smpbellyix;
        depth = depthbelly*10;
        quantileSMP
        medramm(nn) = median(medrammbelly);
        medstrength(nn) = median(medstrengthbelly);
        
        nn = nn + 1;
        % Loop Over VS Files
        smpIx = smpvsix;
        depth = depthvs*10;
        quantileSMP
        medramm(nn) = median(medrammvs);
        medstrength(nn) = median(medstrengthvs);

        plotSondeProfiles
                
        clear ('smprutix','smpbellyix','smpvsix');
    end 
end
%% Median Statistic of RamSonde and SMP for Correlation
% boxRamm = [];
% boxRammIx = [];

boxSMP = [];
boxSMPIx = [];
boxNormSMP = [];
boxSubNormSMP = [];
boxSMP95 = [];
for ii = 1:size(medM,1)
    medSMP(ii,:) = nanmedian(medM{ii}(:,:,1),1);
    SMP5(ii,:) = nanmedian(medM{ii}(:,:,2),1);
    SMP95(ii,:) = nanmedian(medM{ii}(:,:,3),1);
    normMedSMP(ii,:) = nanmedian(normMedM{ii});
%     medramm(ii,1) = median(E(:,rammlocs(ii)));
%     E5(ii,1) = quantile(E(:,rammlocs(ii)),0.05);
%     E95(ii,1) = quantile(E(:,rammlocs(ii)),0.95);
%     boxRamm = [boxRamm;E(:,rammlocs(ii))];
%     boxRammIx = [boxRammIx;rammlocs(ii)*ones(length(E(:,rammlocs(ii))),1)];
    boxSMP = [boxSMP;medM{ii}(:,:,1)];
    boxSMP95 = [boxSMP95;medM{ii}(:,:,3)];
    boxNormSMP = [boxNormSMP;normMedM{ii}];
    boxSubNormSMP = [boxSubNormSMP;subNormMedM{ii}];
    boxSMPIx = [boxSMPIx;ii*ones(size(medM{ii},1),1)];
end
boxSMPgroup = zeros(size(boxSMPIx));
boxrutmem = ismember(boxSMPIx, [1:3:24]);
boxrutIx = find(boxrutmem);
boxbellymem = ismember(boxSMPIx, [2:3:24]);
boxbellyIx = find(boxbellymem);
boxvsmem = ismember(boxSMPIx, [3:3:24]);
boxvsIx = find(boxvsmem);

boxSMPgroup(boxrutIx) = 1;
boxSMPgroup(boxbellyIx) = 2;
boxSMPgroup(boxvsIx) = 3;

%% SMP Box Plot Strength
figure();
        snowLabs = {'Rut','Belly','Virgin Snow'};
        subplot(1,5,[1])
        boxplot(boxSMP(:,2,1),boxSMPgroup,'notch','on','Labels',snowLabs)
%         boxplot(boxSMP95(:,2,1),boxSMPgroup,'notch','on','Labels',snowLabs)
%         boxplot([boxSMP(:,2,1);boxNormSMP(:,2);boxSubNormSMP(:,2)],[boxSMPgroup;3+boxSMPgroup;6+boxSMPgroup],'notch','on')%,'Labels',[snowLabs,snowLabs,snowLabs])
% boxplot(boxNormSMP(:,2),boxSMPgroup,'notch','on','Labels',snowLabs)
%         xlabel('Snow Condition')
        title('Rupture Force')
        xtickangle(-20)
                        set(gca,'fontweight','bold','fontsize',12)
%                                         ylim([0 0.5])

        subplot(1,5,[2])
        q = quantile(boxSMP(:,10,1),[0.25 0.785]);
        q25 = q(1);
        q75 = q(2);
        boxplot(boxSMP(:,10,1),boxSMPgroup,'notch','on','Labels',snowLabs)
%         h = findobj(gcf,'box')
%         set(h(5,1), 'YData', [q25 q75 q75 q25 q75]);% blue box  

%         boxplot(boxSMP95(:,10,1),boxSMPgroup,'notch','on','Labels',snowLabs)
% %         boxplot([boxSMP(:,10,1);boxNormSMP(:,10);boxSubNormSMP(:,10)],[boxSMPgroup;3+boxSMPgroup;6+boxSMPgroup],'notch','on')%,'Labels',[snowLabs,snowLabs,snowLabs])
%         boxplot(boxNormSMP(:,10),boxSMPgroup,'notch','on','Labels',snowLabs)

        title( 'Mean Penetration Force')
%         title( 'L')
%         xlabel('Snow Condition')
                xtickangle(-20)
                set(gca,'fontweight','bold','fontsize',12)
        ylim([-1 10])
        subplot(1,5,[4])
        boxplot(boxSMP(:,7,1),boxSMPgroup,'notch','on','Labels',snowLabs);
%         boxplot(boxSMP95(:,1,1),boxSMPgroup,'notch','on','Labels',snowLabs)
%         boxplot([boxSMP(:,7,1);boxNormSMP(:,7);boxSubNormSMP(:,7)],[boxSMPgroup;3+boxSMPgroup;6+boxSMPgroup],'notch','on')%,'Labels',[snowLabs,snowLabs,snowLabs])
%         boxplot(boxNormSMP(:,7,1),boxSMPgroup,'notch','on','Labels',snowLabs)

%         xlabel('Snow Condition')
                xtickangle(-20)
                set(gca,'fontweight','bold','fontsize',12)
                ylim([0 .5])
        title('Strength')
        subplot(1,5,[3])
        boxplot(boxSMP(:,14,1),boxSMPgroup,'notch','on','Labels',snowLabs)
%         boxplot(boxSMP95(:,14,1),boxSMPgroup,'notch','on','Labels',snowLabs)
%         boxplot([boxSMP(:,14,1);boxNormSMP(:,14);boxSubNormSMP(:,14)],[boxSMPgroup;3+boxSMPgroup;6+boxSMPgroup],'notch','on')%,'Labels',[snowLabs,snowLabs,snowLabs])
%         boxplot(boxNormSMP(:,14,1),boxSMPgroup,'notch','on','Labels',snowLabs)

        title('Density')
        xtickangle(-20)
        ylim([0 500])
        set(gca,'fontweight','bold','fontsize',12)
        subplot(1,5,5)
        boxplot(boxRamm,boxGroup(:,2),'notch','on','Labels',snowLabs)
        title('Ram Hardness')
        ylim([0 150])
                xtickangle(-20)
                set(gca,'fontweight','bold','fontsize',12)
%         ylabel({'Strength'},'rotation',270,'units','normalized','position',[1.25,.5,0])
%         xlabel('Snow Condition')%,'units','normalized','position',[.5,-.25,0])
%         set(gca,'YAxisLocation','right','fontweight','bold','units','normalized','position',[0.75,0.1125,0.155,0.44])
%         ylim([0,1])
        
% % Try All Strengths
% Sruts = cat(1,Sbox{1:3:24});
% Sbelly = cat(1,Sbox{2:3:24});
% Svs = cat(1,Sbox{3:3:24});
% boxS = [Sruts;Sbelly;Svs];
% groupS = [ones(size(Sruts));2*ones(size(Sbelly));3*ones(size(Svs))];
% figure();
%         snowLabs = {'Rut','Belly','Virgin Snow'};
% %         subplot(5,4,[12,20])
%         boxplot(boxS(:),groupS(:),'notch','on','Labels',snowLabs)
% %         xtickangle(-10)
%         ylabel({'Strength'},'rotation',270,'units','normalized','position',[1.25,.5,0])
%         xlabel('Snow Condition')%,'units','normalized','position',[.5,-.25,0])
% %         set(gca,'YAxisLocation','right','fontweight','bold','units','normalized','position',[0.75,0.1125,0.155,0.44])
%         ylim([0,2])

%% Compute Correlation Coefficient
for ii = 1:size(medSMP,2)
    %Use Median for Correlation
    [R,PP,RL,RU] = corrcoef(medSMP(:,ii),medramm);
    % Use Max for Correlation
%     [R,PP,RL,RU] = corrcoef(SMP95(:,ii),medramm);
    % Use Median Normalized for Correlation
%     [R,PP,RL,RU] = corrcoef(normMedSMP(:,ii),medramm);
    XC(1,ii) = R(1,2);
    XC(2,ii) = PP(1,2);
    XC(3,ii) = RL(1,2);
    XC(4,ii) = RU(1,2);
end

% Compute Correlation Coefficient Strength
for ii = 1:size(medSMP,2)
    [R,PP,RL,RU] = corrcoef(medSMP(:,ii),medstrength);
    XCs(1,ii) = R(1,2);
    XCs(2,ii) = PP(1,2);
    XCs(3,ii) = RL(1,2);
    XCs(4,ii) = RU(1,2);
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

    plot(medSMP(:,subarray(ii)),medramm,'ok')
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
    text(0.075,0.775,['R = ', num2str(XC(1,subarray(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.075,0.575,['p = ', num2str(XC(2,subarray(ii)),'%4.2f')],'units','normalized','fontweight','bold')
    
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

%% Plot Strength
% medramm2 = medramm;
legtag = [SVtag{1},'Rut','Belly','Virgin Snow'];
medramm2 = medramm([[1:3:end],[2:3:end],[3:3:end]]);
medSMP2 = medSMP([[1:3:end],[2:3:end],[3:3:end]],:);
symbology = 'oooddd+++***sssxxx^^^ppp';
symbology2 = repmat('od+*sx^p',1,3);
colors = repmat('krb',1,8);
colors2 = 'kkkkkkkkrrrrrrrrbbbbbbbb';
figure();
subz = [1,2,2,2,3,4,5];
for ii = [1,2,5,6,7]
%     figure();
if ii == 7
    subplot(2,3,[5,5.75])
else
    
   subplot(2,3,subz(ii))
end
    h = plot(medSMP(:,ii),medramm,'ok');
    hold on; lsline
    set(h,'visible','off')

    for jj = 1:length(medramm)
        hold on
        %plot(medSMP(jj,ii),medramm(jj),[symbology(jj) colors(jj)],'markersize',10,'markerfacecolor',colors(jj))
        if jj <= 8
        h1(jj) = plot(medSMP2(jj,ii),medramm2(jj),[symbology2(jj) colors2(jj)],'markersize',10,'markerfacecolor',colors2(jj));
        else
            plot(medSMP2(jj,ii),medramm2(jj),[symbology2(jj) colors2(jj)],'markersize',10,'markerfacecolor',colors2(jj));
            ylim([0,150])
        end
        if jj == 8
%             legend(h1,SVtag{1})%,'update','off')
        end
        
    

%     set(gca,'xscale','log')
    end
        h2 = plot(0,0,'k','linewidth',2);
        h3 = plot(0,0,'r','linewidth',2);
        h4 = plot(0,0,'b','linewidth',2);
    set(gca,'fontweight','bold')
        title(titleVar{ii})
%     text(0.675,0.875,['R = ', num2str(CC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.675,0.675,['p = ', num2str(CC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.75,0.2,['R = ', num2str(XC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.75,0.1,['p = ', num2str(XC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
    if ii == 7
        legend([h1,h2,h3,h4],legtag,'position',[[0.85 .1525 0.1 0.2]],'units','normalized')
%         legend([h1,h2,h3,h4],legtag,'location','southeastoutside')

    end
%                 legend([h1,h2,h3,h4],legtag,'location','northwest')%,'update','off')
%                 legend([h2,h3,h4,h1],'Rut','Belly','Virgin Snow',SVtag{1})%,'update','off')

end

% for ii = [1,2,5,6,7]
%     figure();
% % if ii == 7
% %     subplot(2,3,[5,6])
% % else
% %     
% %    subplot(2,3,subz(ii))
% % end
%     h = plot(medSMP(:,ii),medramm,'ok');
%     hold on; lsline
%     set(h,'visible','off')
% 
%     for jj = 1:length(medramm)
%         hold on
%         %plot(medSMP(jj,ii),medramm(jj),[symbology(jj) colors(jj)],'markersize',10,'markerfacecolor',colors(jj))
%         if jj <= 8
%         h1(jj) = plot(medSMP2(jj,ii),medramm2(jj),[symbology2(jj) colors2(jj)],'markersize',10,'markerfacecolor',colors2(jj));
%         else
%             plot(medSMP2(jj,ii),medramm2(jj),[symbology2(jj) colors2(jj)],'markersize',10,'markerfacecolor',colors2(jj));
%             ylim([0,150])
%         end
%         if jj == 8
% %             legend(h1,SVtag{1})%,'update','off')
%         end
%         
%     
% 
% %     set(gca,'xscale','log')
%     end
%         h2 = plot(0,0,'k','linewidth',2);
%         h3 = plot(0,0,'r','linewidth',2);
%         h4 = plot(0,0,'b','linewidth',2);
%     set(gca,'fontweight','bold')
%         title(titleVar{ii})
% %     text(0.675,0.875,['R = ', num2str(CC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
% %     text(0.675,0.675,['p = ', num2str(CC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.875,0.15,['R = ', num2str(XC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.875,0.1,['p = ', num2str(XC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
% %     if ii == 7
% %         legend([h1,h2,h3,h4],legtag,'position',[[0.85 .1525 0.1 0.2]],'units','normalized')
% %         legend([h1,h2,h3,h4],legtag,'location','southeastoutside')
% 
% %     end
%                 legend([h1,h2,h3,h4],legtag,'location','northwest')%,'update','off')
% %                 legend([h2,h3,h4,h1],'Rut','Belly','Virgin Snow',SVtag{1})%,'update','off')
% 
% end
figure();
subz = [1,2,2,2,3,4,5];
for ii = [2,7,14]
%     figure();
if ii == 7
    subplot(3,3,[4,5,5.5])
elseif ii == 2
    
   subplot(3,3,[1,2.5])
else
    subplot(3,3,[7,8.5])
end
    h = plot(medSMP(:,ii),medramm,'ok');
    hold on; lsline
    set(h,'visible','off')

    for jj = 1:length(medramm)
        hold on
        %plot(medSMP(jj,ii),medramm(jj),[symbology(jj) colors(jj)],'markersize',10,'markerfacecolor',colors(jj))
        if jj <= 8
        h1(jj) = plot(medSMP2(jj,ii),medramm2(jj),[symbology2(jj) colors2(jj)],'markersize',10,'markerfacecolor',colors2(jj));
        else
            plot(medSMP2(jj,ii),medramm2(jj),[symbology2(jj) colors2(jj)],'markersize',10,'markerfacecolor',colors2(jj));
            ylim([0,150])
        end
        if jj == 8
%             legend(h1,SVtag{1})%,'update','off')
        end
        
    

%     set(gca,'xscale','log')
    end
        h2 = plot(0,0,'k','linewidth',2);
        h3 = plot(0,0,'r','linewidth',2);
        h4 = plot(0,0,'b','linewidth',2);
    set(gca,'fontweight','bold')
    ylabel('Strength Index')
        title(titleVar{ii})
%     text(0.675,0.875,['R = ', num2str(CC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.675,0.675,['p = ', num2str(CC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.05,0.9,['R = ', num2str(XC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.05,0.8,['p = ', num2str(XC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
    if ii == 7
        legend([h1,h2,h3,h4],legtag,'position',[[0.8375 .1525 0.1 0.2]],'units','normalized')
%         legend([h1,h2,h3,h4],legtag,'location','southeastoutside')

    end
%                 legend([h1,h2,h3,h4],legtag,'location','northwest')%,'update','off')
%                 legend([h2,h3,h4,h1],'Rut','Belly','Virgin Snow',SVtag{1})%,'update','off')

end
%% Figure with 4 5 significant correlations
figure();
subz = [1,2,2,2,3,4,5];
for ii = [2,7,10,14]
%     figure();
if ii == 2
    subplot(2,3,[1])
elseif ii == 7
    
   subplot(2,3,[5])
elseif ii == 10
    subplot(2,3,[2])
elseif ii == 14
    subplot(2,3,4)
% elseif ii == 17
%     subplot(3,2,5)
end
    h = plot(medSMP(:,ii),medramm,'ok');
    hold on; lsline
    set(h,'visible','off')

    for jj = 1:length(medramm)
        hold on
        %plot(medSMP(jj,ii),medramm(jj),[symbology(jj) colors(jj)],'markersize',10,'markerfacecolor',colors(jj))
        if jj <= 8
        h1(jj) = plot(medSMP2(jj,ii),medramm2(jj),[symbology2(jj) colors2(jj)],'markersize',10,'markerfacecolor',colors2(jj));
        else
            plot(medSMP2(jj,ii),medramm2(jj),[symbology2(jj) colors2(jj)],'markersize',10,'markerfacecolor',colors2(jj));
            ylim([0,150])
        end
        if jj == 8
%             legend(h1,SVtag{1})%,'update','off')
        end
        
    

%     set(gca,'xscale','log')
    end
        h2 = plot(0,0,'k','linewidth',2);
        h3 = plot(0,0,'r','linewidth',2);
        h4 = plot(0,0,'b','linewidth',2);
    set(gca,'fontweight','bold')
    ylabel('Hardness Index')
        title(titleVar{ii})
%     text(0.675,0.875,['R = ', num2str(CC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
%     text(0.675,0.675,['p = ', num2str(CC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.05,0.9,['R = ', num2str(XC(1,ii),'%4.2f')],'units','normalized','fontweight','bold')
    text(0.05,0.8,['p = ', num2str(XC(2,ii),'%4.2f')],'units','normalized','fontweight','bold')
    if ii == 7
        legend([h1,h2,h3,h4],legtag,'position',[[0.8375 .1525 0.1 0.2]],'units','normalized')
%         legend([h2,h3,h4],legtag,'position',[[0.8375 .1525 0.1 0.2]],'units','normalized')
%         legend([h1],legtag,'position',[[0.8375 .1525 0.1 0.2]],'units','normalized')

%         columnlegend(3, legtag,'location','SouthEast')
%         legend([h1,h2,h3,h4],legtag,'location','southeastoutside')

    end
%                 legend([h1,h2,h3,h4],legtag,'location','northwest')%,'update','off')
%                 legend([h2,h3,h4,h1],'Rut','Belly','Virgin Snow',SVtag{1})%,'update','off')

end
%% plot profiles


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
%         boxplot(boxramm,boxrammIx,'Labels',SVtag)
%         xtickangle(-10)
%         ylabel('Shear Vane (kPa)','rotation',270,'units','normalized','position',[1.0625,.5,0])
%         xlabel('Location','units','normalized','position',[.5,-.1,0])
%         set(gca,'YAxisLocation','right','fontweight','bold')
% %         suplabel('SMP','x')
% %         suplabel('Shear Vane','y')
%         end
% end
