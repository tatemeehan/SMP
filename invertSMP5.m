function r = invertSMP5(F,z,pfthresh,fthresh,mu,A,theta)
% HPM 04/10/04, 03/04/05, 02/07/07, 05/27/08
% 06/04/08 - now include possibility of a threshhold that is a percentage of
%    max rupture force (pfthresh)
% 06/16/08 - removed fthresh, sthresh - just use pfthresh now, and changed
%    order of input parameters to facilitate default values
% 07/13/08 - include fthresh again, and output Nm,N_T,Na
% this function inverts the SMP signal for mictrostructural parameters 
% NOTE: this calculates Ne and Pc2 using Sturm&al04 method, and gives
%   uncertainty on delta.  If fast calculation is needed, and don't want 
%   Ne,Pc2, or delta uncertainty, use invertSMP
% INPUT: F = [N] penetration force, vertical
%        z = [mm] distance, vertical
% pfthresh = [] fraction of maximum rupture force to threshhold (0-1, dynamic)
% fthresh = [N] minimum rupture force threshold (static) 
%      mu = [] coefficient of friction [0.25]
%        A = [mm^2] cross-sectional area of indenter (pi*r^2) [SMP r=2.5mm]
%   theta = [rad] half-angle of cone.  [SMP=30*pi/180], flat indenter=pi/2 
% OUTPUT: r.L = [mm] structural element length 
%         r.delta = [mm] deflection at rupture normal to tip [Sturm&al,04]
%         r.delta_old = [mm] deflection at rupture normal to tip [J&S,99]
%         r.fn = [N] mean rupture force normal to tip
%         r.N_T = [] total number of ruptures in segment
%         r.Ne = [] number of engaged elements
%         r.Na = [] number of available elements
%         r.Nm = [mm^-1] number of measured ruptures/mm (N_T can be calculated from L or with cal_L2.m)
%         r.N_T = [] number of "true" total ruptures/mm, after mod #1
%         r.Pc = [] Probability of contact = delta/L
%         r.Pc2 = [] Probability of contact = Ne/Na
%         r.k = [N/mm] stiffness = fn/delta
%         r.E = [N/mm^2] elastic modulus = k/L
%         r.sig = [N/mm^2] compressive strength = fn/L^2

% SET DEFAULTS
if nargin<7
    theta=30*pi/180; % default SMP cone half-angle
    if nargin<6
        A=pi*(2.5^2); % default SMP cross-sectional area
        if nargin<5
            mu=0.25; % default coefficent of friction
            if nargin<4
                fthresh=0.014; % Lutz, pers. comm 07/08
                if nargin<3
                    pfthresh=0.1; % minimum force value [10% of max f]
                end
            end
        end
    end
end


% FIRST CORRECT FORCE,DISPLACEMENT FOR GEOMETRY OF TIP
F=F./(sin(theta)+mu*cos(theta)); % Force normal to tip, accounting for friction
z=z*sin(theta); % displacement normal to tip
fthresh=fthresh./(sin(theta)+mu*cos(theta)); % threshold value normal to tip
%%%%%%%%%%%%%%%%%%%%%% NOISE REMOVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find ruptures, threshold ruptures to remove noise, calculate rupture forces
%disp('finding all ruptures...')
[Fall,zall] = get_maxmin2(F,z); % this gets pre-minima, peak val, and post-minima
if length(Fall)==0
    zall=[NaN NaN NaN];
    Fall=[NaN NaN NaN];
end
f_pr=Fall(:,2)-Fall(:,3); % rupture forces
%disp('Noise removal: threshholding ruptures based on rupture slope...')
%f_build=Fall(:,2)-Fall(:,1); % building forces
dF_dz=f_pr./(zall(:,3)-zall(:,2)); % slope of rupture
ind=find(f_pr>(max(f_pr)*pfthresh) & f_pr>fthresh); % dynamic and static threshholding
Fall=Fall(ind,:);  % remove points without this steep slope
zall=zall(ind,:);
f_pr=Fall(:,2)-Fall(:,3); % recalculate rupture forces (same as f_pr=f_pr(ind))
N_T=length(f_pr); % total number of useable peaks 
f_pr = f_pr_correct(N_T,z,zall,F,f_pr); % apply modification #4
if N_T>=1
    % PREALLOCATE
    delta=zeros(N_T-2,1); % displacement at rupture normal to tip
    Ne=zeros(N_T,1); % number of structural elements engaged at rupture
    for i=1:N_T % loop over all peaks, 
        z0=zall(i,2); % reference location
        dd=zall(i:N_T,2)-z0; % delta d values for all subsequent elements
        F_PR=f_pr(i:N_T); % rupture force of current and all subsequent ruptures
        % initialize
        F_Tn=F_PR(1); % current rupture force 
        nse=1; % one element involved
        % iteratively increase number of involved elements until total measured force is reached
        while F_Tn<Fall(i,2)  % while the total estimated force is less than measured total force...
            nse=nse+1; % number of involved grains
            if (nse+1)>length(dd) % if the number of elements involved exceeds the end of the signal
                nse=NaN; % this point will not be calculated
                delta_max=NaN;
                break % stop calculation (break out of while loop)
            else
                ddt=dd(1:nse); % all delta d values (for all involved grains)
                delta_max=dd(nse+1); % largest value delta_nr can take, must be > than last involved element
                F_Tn=sum(F_PR(1:nse).*(1-ddt./delta_max)); 
            end
        end
        % now we know nse, calculate delta
        if ~isfinite(nse) % if we ran out of peaks at end of signal
            delta(i)=NaN;   % make all values NaN
            Ne(i)=NaN;
        elseif nse==1 % if only the current grain active
            delta(i)=zall(i,2)-zall(i,1); % deflection at rupture
            Ne(i)=1;
        else % otherwise delta=delta_max
            % apply modification #2, solve exactly for delta
            num=sum(F_PR(1:nse).*dd(1:nse)); % numerator
            denom=sum(F_PR(1:nse))-Fall(i,2); % denominator
            delta(i)=num/denom; % deflection at rupture - modification #2,3
            Ne(i)=nse; % number of structural elements involved
        end
    end
    delta=delta(find(isfinite(delta))); % get rid of NaNs
    delta=quantile(delta,[0.50]); % displacement at rupture - just use central value
    dz=(max(z)-min(z))/sin(theta); % length of segment, in vertical direction 
    % apply modification #1
    Nm=N_T/dz; % mean number of measured ruptures/mm, rounded to nearest 0.01
    [Ln,N_T]=cal_L2(Nm); % calculate L, and number of true ruptures/mm % modification #1
    N_o=(N_T-Nm)*dz; % number of overlaps in sample
    N_T=N_T*dz; % true total number of ruptures in sample
    fn=sum(f_pr)./N_T; % mean rupture force - depends on N_r
    Ne=Ne(find(isfinite(Ne))); % get rid of NaNs
    Ne=quantile(Ne,[0.5]); % number of elements engaged at rupture - just use central value
    Fm=mean(F); % mean penetration force
    dz2=dz*sin(theta); % length of segment in direction normal to tip
    delta_old=2*dz2*Fm/sum(f_pr); % displacement at rupture, indepenent of Ne ([J&S,99] eq. 22+25)
    r.L=Ln; r.delta=delta; r.delta_old=delta_old; r.fn=fn; r.N_T=N_T; r.Ne=Ne;
    NaL=A/sin(theta)./(r.L.^2); % total number of availble elements
    r.Pc=r.delta./r.L; r.Pc2=r.Ne./NaL; % 
    r.Na=NaL; % total number of available elements
    r.Nm=Nm; % total number of measured ruptures/mm
    r.N_T=N_T; % total number of "true" ruptures/mm after mod #1
    r.k=r.fn./r.delta; % stiffness [N/mm]
    r.E=r.k./r.L; % elastic modulus [N/mm^2]
    r.sig=r.fn./(r.L.^2); % compressive strength [N/mm^2]
else
    disp('no ruptures!!')
    r.L=[NaN NaN NaN]; r.delta=NaN; r.fn=[NaN NaN NaN]; r.N_T=[NaN NaN NaN]; r.Ne=NaN;
    r.Pc=[NaN NaN NaN]; r.Pc2=[NaN NaN NaN]; r.k=[NaN NaN NaN]; r.E=[NaN NaN NaN]; r.sig=[NaN NaN NaN];
    r.Na=[NaN NaN NaN]; r.Nm=NaN;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   % SUBFUNCTIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% get_maxmin2.m
% HPM 03/01/05
% this function gets all maximum, removes duplicate max points, and finds the associated
%   surrounding minima
% INPUT: F = force vector
%        z = depth vector
% OUTPUT: Fall = [Fmin1 Fmax Fmin2] = Force values for pre-min, max, post min
%         zall = [zmin1 zmax zmin2] = associated depths
function [Fall,zall] = get_maxmin2(F,z)
    F=F(:); z=z(:); % force column vectors
    max_index=find(imregionalmax(F)); % this binary function in image processing toolbox finds the indicies to all local max
    % first remove any duplicate max points
    dubs=max_index(find(diff(max_index)==1)); % find two max points that are equal and adjacent (they must be equal to be adjacent)
    Fnew=F; znew=z;
    znew(dubs+1)=mean([znew(dubs), znew(dubs+1)],2); % set the z-value of other point to the midpoint
    Fnew(dubs)=NaN; znew(dubs)=NaN; % set one of duplicate points to NaN
    F=Fnew(~isnan(Fnew)); z=znew(~isnan(znew)); % remove NaN values
    % now lets find all mins and max
    max_index=find(imregionalmax(F)); % this binary function in image processing toolbox finds the indicies to all local max
    min_index=find(imregionalmin(F)); % this bin funct. finds indicies to all local minima
    ind=find(max_index>min(min_index) & max_index<max(min_index)); % make sure there is a minima at the begining and end
    max_index=max_index(ind);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the below method did not work for large data files - meshgrid uses too
    % much memory!!!
    %[Zmin,Zmax]=meshgrid(z(min_index),z(max_index)); % get matricies to compare max/min locations
    %min1_index=sum(Zmin<Zmax,2); % indicies to pre-minima
    %min2_index=min1_index+1; % the next minima will be the post-minima
    %ind=find(min2_index<=length(min_index)); % fix if min1 is last minimia
    %max_index=max_index(ind); min1_index=min1_index(ind); min2_index=min2_index(ind);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    min1_index=zeros(size(max_index)); % preallocate
    min2_index=zeros(size(max_index)); % preallocate
    for n=1:length(max_index)
        %if mod(n,5000)==0
        %   n
        %end
        ind=find(min_index<max_index(n)); % find all minima less than current max
        min1_index(n)=min_index(max(ind)); % index to pre-minima
        min2_index(n)=min_index(max(ind)+1); % index to post minima
    end
    zall=[z(min1_index),z(max_index),z(min2_index)]; % depths of premin,max,postmin
    Fall=[F(min1_index),F(max_index),F(min2_index)]; % penetration force of premin,max,postmin

    %%%%% f_pr_correct - corrects for digitization error during rupture
    
    function f_pr = f_pr_correct(N_T,z,zall,F,f_pr)
        % this function corrects f_pr for increase in F during rupture
        % 05/30/08 - now use dM = actual # measurements during drop, rather
        % than mean
        Nmd=mean((zall(:,3)-zall(:,2))/mean(diff(z))); % mean number of measurements during drop (round up)
        %dz=mean(diff(z)); % depth increment
        for i=1:N_T
            %dM=(zall(i,3)-zall(i,2))/dz; % # measurements during drop
            ind=find(z>zall(i,3)); % find next measurement after min
            ind2=min(ind);
            if length(ind)>3
                if F(ind2+1)>F(ind2) % if next 2 measurements are increasing
                    if F(ind2+2)>F(ind2+1) % if next 3 measurements are increasing
                        f_cor=(F(ind2+2)-F(ind2-1))/3;  % average slope over 3*dz
                    else
                        f_cor=(F(ind2+1)-F(ind2-1))/2; % average slope over 2*dz
                        %   f_cor=F(ind2)-F(ind2-1);
                    end
                else
                    f_cor=F(ind2)-F(ind2-1); % average slope over 1*dz
                end
                f_pr(i)=f_pr(i)+Nmd*f_cor; % apply correction dM times
                %else % if at the end, remove this peak
                %    f_pr(i)=NaN;
            end
        end

     %% cal_L2 -- uses results of MonteCarlo simulation to calculate 
     %%            number of overlaping elements
    function [Ln,N_T]=cal_L2(Nm)
        % this function calculates the structural element length and
        %  true number of ruptures from measured number of ruptures/mm

        A=pi*2.5^2; % base area of cone
        if Nm<62 % changed from 62 to allow calculation...but L values should be flagged below 0.58
%             load NvsL2 % load result from Monte-Carlo (see cal_No10.m)
            Nall = table2array(readtable("NvsL2.csv"));
            L = table2array(readtable("L.csv"));
            ind=find(Nall(:,3)<56); % use only data for Nm<56
            L=L(ind);Nall=Nall(ind,:);
            Ln=[interp1(Nall(:,1),L,Nm,'spline') interp1(Nall(:,3),L,Nm,'spline') interp1(Nall(:,5),L,Nm,'spline')]; % table-lookup for Ln
            N_T=(A./(Ln.^3)); % update true number of ruptures/mm
        else
            %disp([num2str(Nm) ',too many ruptures, L ambiguous!!'])
            Ln=[NaN NaN NaN];
            N_T=[NaN NaN NaN];
        end