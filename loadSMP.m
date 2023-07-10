% loadSMP.m
% HPM 03/24/04 email:hpmarshall@boisestate.edu
% this MATLAB function reads the binary SMP data
%  this was adapted from the IDL code "pntvar400.pro"
% Works only for version 302 SMP data 
% Tate Meehan 07/05/23 - Completed Header Read for Version 509 SMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE:
% >>filename='D:\AVY\SMP_DATA\JOCH031804\040318_Vfeld_SMPRadarNIR\FILE0010.pnt';
% >>d=readSMP(filename)
% Note: this extracts both a force vector and temperature vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: filename - filename,path to SMP binary (.pnt) data file
% OUTPUT: d - structure array containing all data + header info
%         d.filename = given SMP filename
%         d.vers = SMP version number
%         d.comment = information about SMP version
%         d.serial = SMP serial number
%         d.length = total length of the SMP [mm]
%         d.diamter = diamter of the SMP tip [um]
%         d.overload = maximum force load of SMP sensor [N]
%         d.kistlerRange = Force Sensor Range [pC]
%         d.sensitivity = sensor sensitivity [pC/N]
%         d.sensorSerial = Serial Number of Load Cell
%         d.sensorType = binary indication of functionality
%         d.ampRange = amplifier range [pC]
%         d.ampSerial = amplifier Serial Number
%         d.ampType = binary indication of funcitonality for amplifier
%         d.year = year of measurement
%         d.month = month of measurement
%         d.day = day of measurement
%         d.hr = hour of measurement
%         d.min = minute of measurement
%         d.sec = second of measurement
%         d.GPSfix = Quality of GPR Fix
%         d.GPSstate = Statement on GPS Fix identifier
%         d.numsats = number of GPS satelites in view
%         d.pdop = 3D position dilution of precision
%         d.force = [N] vector of force measurements 
%         d.northing = hemisphere of latitude [N or S]
%         d.easting = hemisphere of longitude [E or W]
%         d. latitude = position latitude [ddeg]
%         d.longitude = position longitude [ddeg]
%         d.alitude = position altitude [m]
%         d.batt_V = battery voltage [V]
%         d.vel = average speed [mm/s]
%         d.fsamp = number of force samples
%         d.dzF = distance between force samples [mm]
%         d.zF = depth axis of force data [mm]
%         d.force = raw penetration force [N]
%         d.cF =  conversion factor to force [N/V]
%         d.cP = conversion factor to pressure [MPa/V] 
%         d.tsamp = number of temperature samples
%         d.dzT = distance between temperature samples [mm] 
%         d.zT = depth axis of temperature data [mm]
%         d.temp =  vector of temperature measurements [deg C]
%         d.zero_off = zero-offset to start [mm]
%         d.tempOffset = temperature calibration offset [deg C]
%         d.handOp = binary indication of driver status

function d = loadSMP(filename)
fid = fopen(filename,'r','b');
d.vers = fread(fid,1,'short'); % version
d.nsamp = fread(fid,1,'long'); % number of samples
d.dzF = fread(fid,1,'float'); % distance between force samples [mm]
d.cF = fread(fid,1,'float'); % conversion factor to force (integer to [N])
d.cP = fread(fid,1,'float'); % conversion factor to pressure (integer to [MPa])
d.zero_off = fread(fid,1,'short'); % zero offset to start
d.year = fread(fid,1,'short'); % Time of measurement
d.month = fread(fid,1,'short'); %
d.day = fread(fid,1,'short'); %
d.hr = fread(fid,1,'short'); %
d.min = fread(fid,1,'short'); %
d.sec = fread(fid,1,'short'); %
d.xcoor = fread(fid,1,'double'); % 
d.ycoor = fread(fid,1,'double'); %
d.zcoor = fread(fid,1,'double'); %
d.batt_V = fread(fid,1,'double'); % Battery voltage
d.vel = fread(fid,1,'float'); % average speed [mm/s]
d.loop = fread(fid,1,'long'); % 
d.waypoints = fread(fid,[1 10],'long');% 
d.calstar = fread(fid,[1 10],'short'); % 
d.calend = fread(fid,[1 10],'short'); % 
d.lengthComment = fread(fid,1,'short'); %
d.comment = fread(fid,[1 102],'*char'); %
d.filename = fread(fid,[1 8],'*char'); % smp filename
d.latitude = fread(fid,1,'float'); % [ddeg]
d.longitude = fread(fid,1,'float');% [ddeg]
d.altitude = fread(fid,1,'float')./100;% [m]
d.pdop = fread(fid,1,'float');% precision of dilution [scalar]
d.northing = fread(fid,1,'*char');% North or South Hemisphere
if strcmp(d.northing,'S')
    d.latitude = -d.latitude;
end
d.easting = fread(fid,1,'*char');% East or West Hemisphere
if strcmp(d.easting,'W')
    d.longitude = -d.longitude;
end
[x,y,utmzone] = deg2utm(d.latitude,d.longitude);
d.northing = y;
d.easting = x;
d.numsats = fread(fid,1,'short');% number of GPS satelites in view
d.GPSfix = fread(fid,1,'short');% Quality of GPS fix -1 = NaN, 0 = Not Fixed, 1 = Standalone, 2 = DGPS, 3 = GPS PPS, 4 = RTK fixed, 5 = RTK float
d.GPSstate = fread(fid,1,'*char');%
d.GPSstate = ['Quality of GPS fix -1 = NaN, 0 = Not Fixed, 1 = Standalone, 2 = DGPS, 3 = GPS PPS, 4 = RTK fixed, 5 = RTK float'];
d.reserved1 = fread(fid,1);%
d.xlocal = fread(fid,1,'short');%[ddeg]
d.ylocal = fread(fid,1,'short');%[ddeg]
d.zlocal = fread(fid,1,'short');%[m]
d.thetalocal = fread(fid,1,'short');%[ddeg]
d.reserved2 = fread(fid,[1 62]);%
d.fsamp = fread(fid,1,'long');%
d.tsamp = fread(fid,1,'long');%
d.kistlerRange = fread(fid,1,'short');% [pC]
d.ampRange = fread(fid,1,'short');% [pC]
d.sensitivty = fread(fid,1,'short');% [pC/N]
d.tempOffset = fread(fid,1,"short");% [C]
d.handOp = fread(fid,1,'short');
d.diameter = fread(fid,1,'long');%[um]
d.overload = fread(fid,1,'short');%[N]
d.sensorType = fread(fid,1,'*char');%
d.ampType = fread(fid,1,'*char');%
d.SMPserial = fread(fid,1,'short');% SMP Serial Number
d.length = fread(fid,1,'short');% [mm]
d.reserved3= fread(fid,[1 4]);%
d.sensorSerial = fread(fid,[1 20],'*char');%
d.ampSerial = fread(fid,[1 20],'*char');%
d.reserved4 = fread(fid,[1 80]);%

% Read the Force and Temperature Data
f_block=floor(d.fsamp/256)+1; % number of "256 byte blocks"
t_block=floor(d.tsamp/256)+1; % number of "256 byte blocks"
% now read the data:
status=fseek(fid,512,'bof'); % reposition pointer to start of data
buf=fread(fid,f_block*256+t_block*256,'short'); % total number of data points to read (short integers)
d.force=buf(1:d.fsamp); % force data
temp=buf(f_block*256:f_block*256+d.tsamp); % raw temperature data
temp=temp(1:length(temp)-1);

% now convert the temperature sensor reading -- not making sense, maybe
% sensor is not working.......NOPE, JUST WRONG CALIBRATION CONSTANTS!
a_t= -0.018205007;
b_t= 0.0078710989;
c_t= -0.00098275584;
d_t= 4.2608056e-5;

i=find(temp == 0);
temp(i)=NaN;

logRt=log(temp/0.1/8.0/5.7);
d.temp=1./(a_t + b_t*logRt + c_t*logRt.^2 + d_t*logRt.^3) - 273.15 - 1.9; % temp in 1/100 deg C
d.temp=d.temp/100; % temp in [deg C]
d.force=d.force*d.cF; % force in [N]
d.zF=(0:d.dzF:(d.fsamp-1)*d.dzF)'; % depth scale for force
d.dzT = max(d.zF)./d.tsamp;
d.zT=(0:d.dzT:(d.tsamp-1)*d.dzT)'; % depth scale for force
% d=orderfields(d,[21 3 22 4 1 2 5:20]);
order = [[24,1,23,22,51,52,47,48,42,44,54,49,43,55,50],[[7:12],32,33,31,28:30,25:27,13:15,35:38,19],[16:18],40,2:3,59,57,4,5,41,60,61,58,6,45,46,20,21,34,39,53,56];
d=orderfields(d,order);
d=rmfield(d,{'nsamp','lengthComment','xcoor','ycoor','zcoor','xlocal','ylocal','zlocal','thetalocal','waypoints','loop','calstar','calend','reserved1','reserved2','reserved3','reserved4'});
d.zone = utmzone;
order = [1:26,45,27:44];
d=orderfields(d,order);

status=fclose(fid);
