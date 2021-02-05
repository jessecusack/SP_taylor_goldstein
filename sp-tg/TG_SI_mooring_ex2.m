%% This is the main program to pick out unstable mode families and compute the turbulent energy budget
%
% INPUTS:
% data: a Matlab structure that must contain fileds v, sigma4, z (optional: u, hab)
% u, v : the velocity
% sigma4 : potential density referenced to 4000 dbar
% hab : height above bottom in meters
% z : height in meters
% b : sorted buoyancy profile
% dz: z interval
% D1, D2: only keep data above bottom D1 meters to D2 meters (e.g., D2=1500; D1=0; )
% HOWTO: 1- linear interp; 2-segment mean
% BOT: data profile starts from 1-DEEPEST DATA POINT 2- REAL SEA BOTTOM
% K_init, L_init: initial guesses of wave numbers and wave vector direction
% k_thred: theredhold for k ; 1- 2.*pi./7/(bs/2); 2- do not apply
% bs: bin-size of the measurments (resolution of data)
% tool: 1- the Finite difference method (FD) 2- the Fourier-Galerkin method (FG)
% iBC1, iBCN: boundary condition
%   FD - 
%   (1) velocity: 1=rigid (2nd order: w_{hat}=0; 4th order: w_{hat}=0, w_{hat}_z=0), 0=frictionless (2nd order: w_{hat}_z=kw_{hat}; 4th order: w_{hat}=0, w_{hat}_zz=0)
%   (2) buoyancy: 1=insulating (b_{hat}_z=0 suggesting -\kappa\frac{\pb}{\pz}=0), 0=fixed-buoyancy (b_{hat}=0)                             
%   FG -
%   satisfies impermeable, frictionless, constant-buoyancy boundaries (w_{hat}=0,  w_{hat}_zz=0, b_{hat}=0)
% background viscosity, diffusivity:
%    FD - 
%    nu, kap (constant or (n,1) columns)
%    FG - 
%    Av,Ah,Kv,Kh (constant or (n,1) columns)
%
% OUTPUTS:
% mode families classified by the level of the maximum perturbation kinetic energy
% profiles correspoinding to the direction of mean flow
% 
% CALLS:
% FGM_modefamily.m (change corresponding lines, eg., 83-84 to process the profiles)
% e.g., to sort the buoyancy profile, call DP_1.m; to artifically exclude
% overturns, BBLs, and smooth the data, call DP_2.m; to fill upper layer 
% flow reversals with stagnant and uniform fluid, call DP_3.m: 
%
% FGM_energetics.m (purturbation energy budget)
% FGM_buoyancy.m   (purturbation buoyancy budget)
%
% TO NOTICE:
% code from Smyth: http://salty.oce.orst.edu/wave_analysis/SSF_index.html
% Smyth, W.D., J.N. Moum and J.D. Nash, 2011: ?Narrowband, high-frequency oscillations at the equator. Part II: Properties of shear instabilities", J. Phys. Oceanogr. 41, 412-428. 
% S.Tan, Scripps, 2019/04/15
%
% the example in this script is tested on mooring profiles, 
% only use DP.m and DP_0.m to process data
% use the whole data profile
% FG method (automatically satisfies frictionless BC), real world background viscosity, diffusivity
% use un-sorted sigma4 profiles
% also output energy and buoyancy budgets, k_thred=2
% only support FG method for this script
% adapt the idea of seeking the mode families from Smyth, W. D., Moum, J. N., Li, L., & Thorpe, S. A. (2013). Diurnal Shear Instability, the Descent of the Surface Shear Layer, and the Deep Cycle of Equatorial Turbulence, Journal of Physical Oceanography, 43(11), 2432-2455
% dz = 4 will take about 20 mins to run for a single profile, approximately
    % 3 days for T6 mooring... dz = 5 will significantly reduce the time it
    % takes, dz = 4 centrainly have less errors but the results from dz = 5 is
    % not too bad...
% S.Tan, IOCAS, 2021/02/05

clear all
close all
clc

%% datapath & toolpath & etc. 
% main directory (change to your directory)
mdirec='/Users/tantanmeow/GitHub/SP_taylor_goldstein/sp-tg/';
% data directory (change to the data directory)
datadirec='/Users/tantanmeow/Desktop/WORK/2018-2019/Jesse/proc_data/';
% add Matlab tools to path
addpath(genpath(strcat(mdirec,'tools/')));
% add codes tools to path
addpath(genpath(strcat(mdirec,'codes/')));

%% parameters for the instability scan
tool = 2; % 1- FD ; 2- FG
% ----------  used in SSF.m or vTG.m  ----------  
% rangeing wave number for the mode scan
k=[-1.8:.02:-.6]; K=10.^(k);  %successive values differing by a factor 10.^(.02)=1.0471
L = [0:10:180];%theta
disp(strcat('wave length ranges from ', num2str(2*pi/K(end)), 'm', 'to', num2str(2*pi/10.^(k(1)))))
% wave length ranges from25.0138mto396.4422

% (FD)
% boundary condition
% (1) velocity: 1=rigid (2nd order: w_{hat}=0; 4th order: w_{hat}=0, w_{hat}_z=0), 0=frictionless (2nd order: w_{hat}_z=kw_{hat}; 4th order: w_{hat}=0, w_{hat}_zz=0)
% (2) buoyancy: 1=insulating (b_{hat}_z=0 suggesting -\kappa\frac{\pb}{\pz}=0), 0=fixed-buoyancy (b_{hat}=0)                             
iBC1 = [0,0]; % at z=z(0)     
iBCN = iBC1;  % at z=z(N+1)

% (FG) no need to specify BC : 
% Lian et al. (2020): This choice of basis functions is natural for impermeable, frictionless, constant-buoyancy
% boundaries, and may be adapted for insulating boundaries by expanding b_{hat} in terms of cosines
% rather than sines. But the method is not well-suited for rigid boundaries.

% ----------  used in DP.m  ----------  
% processed data grids (see HOWTO and BOT)
% (to note: a finer grid may better capture the critical level and you may want to use the new Matlab codes to do that since the old codes can be time consuming)
dz=4;
HOWTO = 1;% 1- linear interp; 2-segment mean
BOT=1; % data profile starts from 1-DEEPEST DATA POINT 2- REAL SEA BOTTOM

% ----------  used in FGM_modefamily.m  ----------  
% bin-size of the measurments (resolution of data)
bs=8; % 10m or 16m? I'm not sure

% theredhold for k
k_thred = 2; % 1- 2.*pi./7/(bs/2); 2- do not apply

% butterworth low-pass filter window
fl_bw = 20; % this is because the thrope scale is pretty much this level

% experiment
expri = 2;

% GR_lm: ignore growth rates that are shorter than this (hr^-1)
GR_lm = 1;

% imode: how many modes to retain for the first scan
imode = 10;

% how many critical levels to retain and at one critical level how many
% families to retain
cl_nu = 50;
fa_nu = 50;

newcolors = {'#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k','#F00','#F80','m','#0B0','#00F','#50F','#A0F','#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F','k'};

%% below is an example for the sorted mooring profiles at each time step
% load data
filename=[dir(strcat(datadirec,'T2*')) dir(strcat(datadirec,'T6*'))];
for filenum=1:length(filename)
    load(strcat(datadirec,filename(filenum).name))
    s = datenum; %number of profiles
    % Convert density to buoyancy
    b=nan(size(sig4));
    g=9.81;
    for i = 1:length(s)
        l=find(~isnan(sig4(:,i)));
        b(l,i)=-g*(sig4(l,i)-nanmean(sig4(:,i),1))/nanmean(sig4(:,i),1);
    end
    
    buoyancy = b;
    vdata = v; udata = u; zdata = z; 
    eps_=eps;

    GR=nan(length(s), cl_nu, fa_nu);  % growth rate
    CR=nan(length(s), cl_nu, fa_nu);  % real phase speed
    CI=nan(length(s), cl_nu, fa_nu);  % imaginary phase speed
    CL=nan(length(s), cl_nu, fa_nu); % critical level
    ERR_K=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    Kt=nan(length(s), cl_nu, fa_nu);  % wavenumber
    Phi=nan(length(s), cl_nu, fa_nu);  % wavevector
    % properties of perturbation kinetic energy budget
    Knetm=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    SPm=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    BFm=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    EFnm=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    cEFnm=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    epsm=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    RHm=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    LHm=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy

    Knet=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    SP=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    BF=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    EFn=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    cEFn=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    eps=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    RH=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy
    LH=nan(length(s), cl_nu, fa_nu); % error in perturbation kinetic energy

    END=nan(length(s), cl_nu, fa_nu); % whether the FGM is between two modes with smaller gr but larger and smaller wavelengths

    % flow properties at the critical level
    Vz=nan(length(s), cl_nu, fa_nu);  
    N2=nan(length(s), cl_nu, fa_nu);  
    Ri=nan(length(s), cl_nu, fa_nu);  

    % flow properties of the mean flow
    botz=nan(length(s),1);% sea bottom depth
    hab=nan(length(s),1);% last measurment above bottom
    zw=[-6000:dz:0]';
    % velocity profiles correspond to the angle of the mean flow
    Vm=nan(length(zw),length(s));Vzm=nan(length(zw),length(s));Vzzm=nan(length(zw),length(s));
    % buoyancy, N squre and Richardson number profiles
    Bm=nan(length(zw),length(s));N2m=nan(length(zw),length(s));Rim=nan(length(zw),length(s));
    Phim = nan(length(s),1); % angle of mean flow
    
    NOS = nan(length(s),2); 
for profile=1:length(s)
    data.v=vdata(find(~isnan(vdata(:,profile))),profile);data.u=udata(find(~isnan(vdata(:,profile))),profile);data.z=zdata(find(~isnan(vdata(:,profile))),profile);
    data.sigma4=sig4(find(~isnan(vdata(:,profile))),profile); %use UN-SORTED profiles 
    data.b=buoyancy(find(~isnan(vdata(:,profile))),profile); 
    if ~isempty(data.v) && data.z(1)-data.z(end)>10*dz
        data.hab=HAB(find(zdata(:,profile)==data.z(end)),profile);
        hab(profile)=data.hab;
        if isfield(data, 'hab')
            botz(profile)=data.z(end)-data.hab;
        end
        % D2=[];D1=[];% only keep data above bottom D1 meters to D2 meters
        D2=data.z(1)-botz(profile);D1=0; %whole data profile

        % low-pass filter
        if fl_bw>0
            DZ=diff(data.z);
            DZ=abs(mean(DZ));
            data.v(~isnan(data.v))=butterworth_lp(data.v(~isnan(data.v)),fl_bw/DZ,4,length(data.v(~isnan(data.v))));
            data.u(~isnan(data.u))=butterworth_lp(data.u(~isnan(data.u)),fl_bw/DZ,4,length(data.u(~isnan(data.u))));
            data.b(~isnan(data.b))=butterworth_lp(data.b(~isnan(data.b)),fl_bw/DZ,4,length(data.b(~isnan(data.b))));
        end      

        [v,vz,vzz,b,n2,Ri_r,zz]=DP(data,dz,D1,D2,L(1),HOWTO,BOT);
        [v,vz,vzz,b,n2,Ri_r,zz]=DP_0(v,vz,vzz,b,n2,Ri_r,zz);
        
        % compute kappa following Osborn
        kappa = 0.2.*eps_(:,profile)./N2_overturn(:,profile);
        Kv = interp1(z(:,profile),kappa,zz);
        Kv(find(Kv<1e-6)) = nan;
        Kv(isnan(Kv)) = 1e-6; % background kappa
        Av = Kv; % Pr = 1
        Ah=Av;Kh=Kv;
        kap=Kv;nu=Av;

        % scan the (k,l) plane to pick out FGM families
        figure(profile)
        [GR_family,CR_family,CI_family,CL_family,ERR_family,K_family,L_family,Knetm_family,SPm_family,BFm_family,EFnm_family,cEFnm_family,epsm_family,RHm_family,LHm_family,Knet_family,SP_family,BF_family,EFn_family,cEFn_family,eps_family,RH_family,LH_family,ends_family]=FGM_modefamily(data,dz,D1,D2,HOWTO,BOT,K,L,Av,Ah,Kv,Kh,GR_lm,imode);
        if ~isempty(GR_family)
            print('-djpeg',[mdirec strcat('results/', filename(filenum).name(1:end-4), '_ex',num2str(expri),'_1_',num2str(profile))])
            close

            % sort by growth rate
            GR_family(isnan(GR_family)) = -999;
            GR_family(GR_family==0) = -999;
            [GR_family, index_family] = sort(GR_family,2,'descend'); 
            for i=1:length(GR_family(:,1))
                CR_family(i,:) = CR_family(i,index_family(i,:));
                CI_family(i,:) = CI_family(i,index_family(i,:));
                CL_family(i,:) = CL_family(i,index_family(i,:));
                ERR_family(i,:) = ERR_family(i,index_family(i,:));
                K_family(i,:) = K_family(i,index_family(i,:));
                L_family(i,:) = L_family(i,index_family(i,:));
                Knet_family(i,:) = Knet_family(i,index_family(i,:));
                SP_family(i,:) = SP_family(i,index_family(i,:));
                BF_family(i,:) = BF_family(i,index_family(i,:));
                EFn_family(i,:) = EFn_family(i,index_family(i,:));
                cEFn_family(i,:) = cEFn_family(i,index_family(i,:));
                eps_family(i,:) = eps_family(i,index_family(i,:));
                RH_family(i,:) = RH_family(i,index_family(i,:));
                LH_family(i,:) = LH_family(i,index_family(i,:));
                Knetm_family(i,:) = Knetm_family(i,index_family(i,:));
                SPm_family(i,:) = SPm_family(i,index_family(i,:));
                BFm_family(i,:) = BFm_family(i,index_family(i,:));
                EFnm_family(i,:) = EFnm_family(i,index_family(i,:));
                cEFnm_family(i,:) = cEFnm_family(i,index_family(i,:));
                epsm_family(i,:) = epsm_family(i,index_family(i,:));
                RHm_family(i,:) = RHm_family(i,index_family(i,:));
                LHm_family(i,:) = LHm_family(i,index_family(i,:));
                ends_family(i,:) = ends_family(i,index_family(i,:));
            end
            [a, index_family] = sort(GR_family(:,1),'descend'); 
            GR_family = GR_family(index_family,:);
            CR_family = CR_family(index_family,:);
            CI_family = CI_family(index_family,:);
            CL_family = CL_family(index_family,:);
            ERR_family = ERR_family(index_family,:);
            K_family = K_family(index_family,:);
            L_family = L_family(index_family,:);
            Knet_family = Knet_family(index_family,:);
            SP_family = SP_family(index_family,:);
            BF_family = BF_family(index_family,:);
            EFn_family = EFn_family(index_family,:);
            cEFn_family = cEFn_family(index_family,:);
            eps_family = eps_family(index_family,:);
            RH_family = RH_family(index_family,:);
            LH_family = LH_family(index_family,:);
            Knetm_family = Knetm_family(index_family,:);
            SPm_family = SPm_family(index_family,:);
            BFm_family = BFm_family(index_family,:);
            EFnm_family = EFnm_family(index_family,:);
            cEFnm_family = cEFnm_family(index_family,:);
            epsm_family = epsm_family(index_family,:);
            RHm_family = RHm_family(index_family,:);
            LHm_family = LHm_family(index_family,:);
            ends_family = ends_family(index_family,:);
            GR_family(find(GR_family<0)) = nan;

            % compute shear, N2, and Ri at the level where extract energy (our defined "cl"s)
            Vz_family = nan(size(GR_family));
            N2_family = nan(size(GR_family));
            Ri_family = nan(size(GR_family));
            figure(profile+50)
            set(gcf, 'PaperPosition',[0 0 50 10]);    
            for i=1:length(GR_family(:,1))
                subplot(1,length(GR_family(:,1))+1,i)
                for j=1:length(GR_family(1,:))
                    if ~isnan(GR_family(i,j))
                        angle = L_family(i,j);
                        [v,vz,vzz,b,n2,Ri_r,zz]=DP(data,dz,D1,D2,angle,HOWTO,BOT);
                        [v,vz,vzz,b,n2,Ri_r,zz]=DP_0(v,vz,vzz,b,n2,Ri_r,zz);
                        hold on
                        plot(n2,zz,'b')
                        hold on
                        plot(vz.^2./4,zz,'r')
                        hold on
                        plot([nanmin([vz.^2./4; n2]) nanmax([vz.^2./4; n2])],[CL_family(i,j) CL_family(i,j)],'k','linewidth',1.5)
                        grid on
                        Vz_family(i,j) = interp1(zz,vz,CL_family(i,j));
                        N2_family(i,j) = interp1(zz,n2,CL_family(i,j));
                        Ri_family(i,j) = interp1(zz,Ri_r,CL_family(i,j));
                    end
                end
                title(strcat('max(gr)=',num2str(round(nanmax(GR_family(i,:))*3600*100)/100)))
            end
            subplot(1,length(GR_family(:,1))+1,length(GR_family(:,1))+1)
            plot(log10(Av),zz,'linewidth',1.5);grid on
            title('Av=Kv')
            print('-djpeg',[mdirec strcat('results/', filename(filenum).name(1:end-4), '_ex',num2str(expri),'_2_',num2str(profile))])
            close

            % compute flow properties in the direction of the mean flow
            u_m = trapz(data.z,data.u)/sum(diff(data.z));
            v_m = trapz(data.z,data.v)/sum(diff(data.z));
            Phim(profile) = 180*atan2(v_m,u_m)/pi;
            [v,vz,vzz,b,n2,Ri_r,zz]=DP(data,dz,D1,D2,Phim(profile),HOWTO,BOT);
            [v,vz,vzz,b,n2,Ri_r,zz]=DP_0(v,vz,vzz,b,n2,Ri_r,zz);
            Vm(:,profile)=interp1(zz,v,zw);Vzm(:,profile)=interp1(zz,vz,zw);Vzzm(:,profile)=interp1(zz,vzz,zw);
            Bm(:,profile)=interp1(zz,b,zw);N2m(:,profile)=interp1(zz,n2,zw);
            Rim(:,profile)=interp1(zz,Ri_r,zw);

        %     figure(profile+100)
        %     DP_pic(Vm(:,profile),Vzm(:,profile),Vzzm(:,profile),Bm(:,profile),N2m(:,profile),Rim(:,profile),zw)
        %     set(gcf, 'PaperPosition',[0 0 8 4]);
        %     print('-djpeg',[mdirec strcat('codes_theta/', filename(filenum).name(1:end-4), '_ex',num2str(expri),'_4_',num2str(profile))])

        %    disp(strcat('No of critical levels: ----- ',num2str(length(GR_family(:,1))),' -----'))
        %    disp(strcat('No of modes: ----- ',num2str(length(GR_family(1,:))),' -----'))
            NOS(profile,1)=length(GR_family(:,1));
            NOS(profile,2)=length(GR_family(1,:));

            % plot growth rates
            figure(profile+150)
            set(gcf, 'PaperPosition',[0 0 40 16]);    
            subplot(241)
            for i=1:length(GR_family(:,1))
                hold on
                plot(GR_family(i,:).*3600,CL_family(i,:),'o','linewidth',1.5,'markersize',5)
                xlabel('growth rate (/hr)')
                ylabel('critical level (m)')
                grid on
            end
            subplot(243)
            for i=1:length(GR_family(:,1))
                hold on
                plot(2.*pi./K_family(i,:),L_family(i,:),'o','linewidth',1.5,'markersize',5)
                xlabel('wavelength (m)')
                ylabel('direction (degree)')
                grid on
            end
            hold on
            plot(2.*pi./[nanmin(nanmin(K_family)),nanmax(nanmax(K_family))],[Phim(profile),Phim(profile)],'color',[.5 .5 .5],'linewidth',1.5)
            ylim([0,180])
            set(gca,'ytick',[0,90,180])
            subplot(242)
            for i=1:length(GR_family(:,1))
                hold on
                plot(2.*pi./K_family(i,:),CL_family(i,:),'o','linewidth',1.5,'markersize',5)
                xlabel('wavelength (m)')
                ylabel('critical level (m)')
                grid on
            end
            subplot(246)
            for i=1:length(GR_family(:,1))
                hold on
                plot(L_family(i,:),CL_family(i,:),'o','linewidth',1.5,'markersize',5)
                xlabel('direction (degree)')
                ylabel('critical level (m)')
                grid on
            end
            hold on
            plot([Phim(profile),Phim(profile)],[nanmin(nanmin(CL_family)),nanmax(nanmax(CL_family))],'color',[.5 .5 .5],'linewidth',1.5)
            xlim([0,180])
            set(gca,'xtick',[0,90,180])
            subplot(245)
            for i=1:length(GR_family(:,1))
                hold on
                plot(CR_family(i,:),CL_family(i,:),'o','linewidth',1.5,'markersize',5)
                xlabel('phase speed (m/s)')
                ylabel('critical level (m)')
                grid on
            end
            subplot(247)
            for i=1:length(GR_family(:,1))
                hold on
                plot(ERR_family(i,:),CL_family(i,:),'o','linewidth',1.5,'markersize',5)
                xlabel('KE err')
                ylabel('critical level (m)')
                grid on
            end
            subplot(244)
            for i=1:length(GR_family(:,1))
                hold on
                plot(SPm_family(i,:),CL_family(i,:),'o','linewidth',1.5,'markersize',5)
                xlabel('shear production')
                ylabel('critical level (m)')
                grid on
            end
            subplot(248)
            for i=1:length(GR_family(:,1))
                hold on
                plot(BFm_family(i,:),CL_family(i,:),'o','linewidth',1.5,'markersize',5)
                xlabel('buoyancy flux')
                ylabel('critical level (m)')
                grid on
            end
            print('-djpeg',[mdirec strcat('results/', filename(filenum).name(1:end-4), '_ex',num2str(expri),'_3_',num2str(profile))])
            close

            if length(GR_family(:,1))<=cl_nu&&length(GR_family(1,:))<=fa_nu
                GR(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=GR_family;
                CR(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=CR_family;
                CI(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=CI_family;
                CL(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=CL_family;
                ERR_K(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=ERR_family;
                Kt(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=K_family;
                Phi(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=L_family;
                Knet(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=Knet_family;
                SP(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=SP_family;
                BF(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=BF_family;
                EFn(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=EFn_family;
                cEFn(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=cEFn_family;
                eps(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=eps_family;
                RH(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=RH_family;
                LH(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=LH_family;
                Knetm(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=Knetm_family;
                SPm(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=SPm_family;
                BFm(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=BFm_family;
                EFnm(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=EFnm_family;
                cEFnm(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=cEFnm_family;
                epsm(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=epsm_family;
                RHm(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=RHm_family;
                LHm(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=LHm_family;
                END(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=ends_family;
                Vz(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=Vz_family;
                N2(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=N2_family;
                Ri(profile,1:length(GR_family(:,1)),1:length(GR_family(1,:)))=Ri_family;
            elseif length(GR_family(:,1))>cl_nu&&length(GR_family(1,:))<=fa_nu 
                GR(profile,:,1:length(GR_family(1,:)))=GR_family(1:cl_nu,:);
                CR(profile,:,1:length(GR_family(1,:)))=CR_family(1:cl_nu,:);
                CI(profile,:,1:length(GR_family(1,:)))=CI_family(1:cl_nu,:);
                CL(profile,:,1:length(GR_family(1,:)))=CL_family(1:cl_nu,:);
                ERR_K(profile,:,1:length(GR_family(1,:)))=ERR_family(1:cl_nu,:);
                Kt(profile,:,1:length(GR_family(1,:)))=K_family(1:cl_nu,:);
                Phi(profile,:,1:length(GR_family(1,:)))=L_family(1:cl_nu,:);
                Knet(profile,:,1:length(GR_family(1,:)))=Knet_family(1:cl_nu,:);
                SP(profile,:,1:length(GR_family(1,:)))=SP_family(1:cl_nu,:);
                BF(profile,:,1:length(GR_family(1,:)))=BF_family(1:cl_nu,:);
                EFn(profile,:,1:length(GR_family(1,:)))=EFn_family(1:cl_nu,:);
                cEFn(profile,:,1:length(GR_family(1,:)))=cEFn_family(1:cl_nu,:);
                eps(profile,:,1:length(GR_family(1,:)))=eps_family(1:cl_nu,:);
                RH(profile,:,1:length(GR_family(1,:)))=RH_family(1:cl_nu,:);
                LH(profile,:,1:length(GR_family(1,:)))=LH_family(1:cl_nu,:);
                Knetm(profile,:,1:length(GR_family(1,:)))=Knetm_family(1:cl_nu,:);
                SPm(profile,:,1:length(GR_family(1,:)))=SPm_family(1:cl_nu,:);
                BFm(profile,:,1:length(GR_family(1,:)))=BFm_family(1:cl_nu,:);
                EFnm(profile,:,1:length(GR_family(1,:)))=EFnm_family(1:cl_nu,:);
                cEFnm(profile,:,1:length(GR_family(1,:)))=cEFnm_family(1:cl_nu,:);
                epsm(profile,:,1:length(GR_family(1,:)))=epsm_family(1:cl_nu,:);
                RHm(profile,:,1:length(GR_family(1,:)))=RHm_family(1:cl_nu,:);
                LHm(profile,:,1:length(GR_family(1,:)))=LHm_family(1:cl_nu,:);
                END(profile,:,1:length(GR_family(1,:)))=ends_family(1:cl_nu,:);
                Vz(profile,:,1:length(GR_family(1,:)))=Vz_family(1:cl_nu,:);
                N2(profile,:,1:length(GR_family(1,:)))=N2_family(1:cl_nu,:);
                Ri(profile,:,1:length(GR_family(1,:)))=Ri_family(1:cl_nu,:);
            elseif length(GR_family(:,1))<=cl_nu&&length(GR_family(1,:))>fa_nu 
                GR(profile,1:length(GR_family(:,1)),:)=GR_family(:,1:fa_nu);
                CR(profile,1:length(GR_family(:,1)),:)=CR_family(:,1:fa_nu);
                CI(profile,1:length(GR_family(:,1)),:)=CI_family(:,1:fa_nu);
                CL(profile,1:length(GR_family(:,1)),:)=CL_family(:,1:fa_nu);
                ERR_K(profile,1:length(GR_family(:,1)),:)=ERR_family(:,1:fa_nu);
                Kt(profile,1:length(GR_family(:,1)),:)=K_family(:,1:fa_nu);
                Phi(profile,1:length(GR_family(:,1)),:)=L_family(:,1:fa_nu);
                Knet(profile,1:length(GR_family(:,1)),:)=Knet_family(:,1:fa_nu);
                SP(profile,1:length(GR_family(:,1)),:)=SP_family(:,1:fa_nu);
                BF(profile,1:length(GR_family(:,1)),:)=BF_family(:,1:fa_nu);
                EFn(profile,1:length(GR_family(:,1)),:)=EFn_family(:,1:fa_nu);
                cEFn(profile,1:length(GR_family(:,1)),:)=cEFn_family(:,1:fa_nu);
                eps(profile,1:length(GR_family(:,1)),:)=eps_family(:,1:fa_nu);
                RH(profile,1:length(GR_family(:,1)),:)=RH_family(:,1:fa_nu);
                LH(profile,1:length(GR_family(:,1)),:)=LH_family(:,1:fa_nu);
                Knetm(profile,1:length(GR_family(:,1)),:)=Knetm_family(:,1:fa_nu);
                SPm(profile,1:length(GR_family(:,1)),:)=SPm_family(:,1:fa_nu);
                BFm(profile,1:length(GR_family(:,1)),:)=BFm_family(:,1:fa_nu);
                EFnm(profile,1:length(GR_family(:,1)),:)=EFnm_family(:,1:fa_nu);
                cEFnm(profile,1:length(GR_family(:,1)),:)=cEFnm_family(:,1:fa_nu);
                epsm(profile,1:length(GR_family(:,1)),:)=epsm_family(:,1:fa_nu);
                RHm(profile,1:length(GR_family(:,1)),:)=RHm_family(:,1:fa_nu);
                LHm(profile,1:length(GR_family(:,1)),:)=LHm_family(:,1:fa_nu);
                END(profile,1:length(GR_family(:,1)),:)=ends_family(:,1:fa_nu);
                Vz(profile,1:length(GR_family(:,1)),:)=Vz_family(:,1:fa_nu);
                N2(profile,1:length(GR_family(:,1)),:)=N2_family(:,1:fa_nu);
                Ri(profile,1:length(GR_family(:,1)),:)=Ri_family(:,1:fa_nu);
            elseif length(GR_family(:,1))>cl_nu&&length(GR_family(1,:))>fa_nu 
                GR(profile,:,:)=GR_family(1:cl_nu,1:fa_nu);
                CR(profile,:,:)=CR_family(1:cl_nu,1:fa_nu);
                CI(profile,:,:)=CI_family(1:cl_nu,1:fa_nu);
                CL(profile,:,:)=CL_family(1:cl_nu,1:fa_nu);
                ERR_K(profile,:,:)=ERR_family(1:cl_nu,1:fa_nu);
                Kt(profile,:,:)=K_family(1:cl_nu,1:fa_nu);
                Phi(profile,:,:)=L_family(1:cl_nu,1:fa_nu);
                Knet(profile,:,:)=Knet_family(1:cl_nu,1:fa_nu);
                SP(profile,:,:)=SP_family(1:cl_nu,1:fa_nu);
                BF(profile,:,:)=BF_family(1:cl_nu,1:fa_nu);
                EFn(profile,:,:)=EFn_family(1:cl_nu,1:fa_nu);
                cEFn(profile,:,:)=cEFn_family(1:cl_nu,1:fa_nu);
                eps(profile,:,:)=eps_family(1:cl_nu,1:fa_nu);
                RH(profile,:,:)=RH_family(1:cl_nu,1:fa_nu);
                LH(profile,:,:)=LH_family(1:cl_nu,1:fa_nu);
                Knetm(profile,:,:)=Knetm_family(1:cl_nu,1:fa_nu);
                SPm(profile,:,:)=SPm_family(1:cl_nu,1:fa_nu);
                BFm(profile,:,:)=BFm_family(1:cl_nu,1:fa_nu);
                EFnm(profile,:,:)=EFnm_family(1:cl_nu,1:fa_nu);
                cEFnm(profile,:,:)=cEFnm_family(1:cl_nu,1:fa_nu);
                epsm(profile,:,:)=epsm_family(1:cl_nu,1:fa_nu);
                RHm(profile,:,:)=RHm_family(1:cl_nu,1:fa_nu);
                LHm(profile,:,:)=LHm_family(1:cl_nu,1:fa_nu);
                END(profile,:,:)=ends_family(1:cl_nu,1:fa_nu);
                Vz(profile,:,:)=Vz_family(1:cl_nu,1:fa_nu);
                N2(profile,:,:)=N2_family(1:cl_nu,1:fa_nu);
                Ri(profile,:,:)=Ri_family(1:cl_nu,1:fa_nu);
            end
        end
    end
end
date = datenum;
eval(strcat('save TG_SI_',  filename(filenum).name(1:end-4), '_ex', num2str(expri), '.mat date lon lat botz hab zw GR CR CI CL ERR_K K L Kt Phi Knetm SPm BFm EFnm cEFnm epsm RHm LHm Knet SP BF EFn cEFn eps RH LH END Vz N2 Ri Vm Vzm Vzzm Bm N2m Rim Phim NOS'))
end

% sp = reshape(squeeze(SP(profile,:,:)),50*50,1);
% bf = reshape(squeeze(BF(profile,:,:)),50*50,1);
% gr = reshape(squeeze(GR(profile,:,:)),50*50,1);
% err = reshape(squeeze(ERR_K(profile,:,:)),50*50,1);
% cl = reshape(squeeze(CL(profile,:,:)),50*50,1);
% 
% figure
% subplot(221)
% histogram(gr(~isnan(gr)),10);title('gr')
% subplot(222)
% histogram(sp(~isnan(gr)),10);title('sp (blue), bf (red)')
% hold on
% histogram(bf(~isnan(gr)),10);
% subplot(223)
% histogram(err(~isnan(gr)),10);title('err')
% subplot(224)
% histogram(cl(~isnan(gr)),10);title('cl')
% 
% disp('binsize = 5 m')
% disp(strcat('number of modes = ', num2str(length(gr(~isnan(gr))))))
% disp(strcat('number of SI modes = ', num2str(length(find((sp(~isnan(gr))>0&bf(~isnan(gr))<0) | (sp(~isnan(gr))>0&bf(~isnan(gr))>0&sp(~isnan(gr))>bf(~isnan(gr))))))))
% disp(strcat('number of CI modes = ', num2str(length(find((sp(~isnan(gr))<0&bf(~isnan(gr))>0) | (sp(~isnan(gr))>0&bf(~isnan(gr))>0&sp(~isnan(gr))<bf(~isnan(gr))))))))

% binsize = 4 m
% number of modes =103
% number of SI modes =98
% number of CI modes =5
% 
% binsize = 5 m
% number of modes =105
% number of SI modes =92
% number of CI modes =13

% binsize = 10 m
% number of modes =42
% number of SI modes =38
% number of CI modes =4