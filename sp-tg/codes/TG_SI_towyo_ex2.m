%% This is the main program to compute fastest growth rate of shear instabilities
%
% INPUTS:
% data: a Matlab structure that must contain fileds v, sigma4, z (optional: u, hab)
% u, v : the velocity
% sigma4 : potential density referenced to 4000 dbar
% hab : height above bottom in meters
% z : height in meters
% b : sorted buoyancy profile
% dz: z interval
% D1, D2: only keep data above bottom D1 meters to D2 meters (e.g., D2=1500 or interface height;D1=0; )
% HOWTO: 1- linear interp; 2-segment mean
% BOT: data profile starts from 1-DEEPEST DATA POINT 2- REAL SEA BOTTOM
% K_init, L_init: initial guesses of wave numbers in x and y directions
% k_thred: theredhold for k ; 1- 2.*pi./7/(bs/2); 2- do not apply
% bs: bin-size of the measurments (resolution of data)
% tool: 1- the Finite difference method (FD) 2- the Fourier-Galerkin method (FG)
% iBC1, iBCN: boundary condition
%   FD - 
%   (1) velocity: 1=rigid (2nd order: w_{hat}=0; 4th order: w_{hat}=0, w_{hat}_z=0), 0=frictionless (2nd order: w_{hat}_z=kw_{hat}; 4th order: w_{hat}=0, w_{hat}_zz=0)
%   (2) buoyancy: 1=insulating (b_{hat}_z=0 suggesting -\kappa\frac{\pb}{\pz}=0), 0=fixed-buoyancy (b_{hat}=0)                             
%   FG -
%   satisfies impermeable, frictionless, constant-buoyancy boundaries (w_{hat}=0,  w_{hat}_zz=0, b_{hat}=0)
% background viscosity, diffusivity
%    FD - 
%    nu, kap (constant or (n,1) columns)
%    FG - 
%    Av,Ah,Kv,Kh (constant or (n,1) columns)
% scan2: 1- only scan once; 2- perform 2nd scan
%
% OUTPUTS:
% KFGM & GR - wave number & fastest growth rate
%
% CALLS:
% FGM_towyo_ex1.m (change corresponding lines, eg., 60-68 to process the profiles)
% e.g., to sort the buoyancy profile, call DP_1.m; to artifically exclude
% overturns, BBLs, and smooth the data, call DP_2.m; to fill upper layer 
% flow reversals with stagnant and uniform fluid, call DP_3.m: 
%
% SI_pic_1(sigs,I,K,i,cphs,inst);
% SI_pic_2(w,K,i,zz);
% DP_pic_1(v,vz,vzz,b,n2,Ri_r,zz,KFGM(profile));
%
% TO NOTICE:
% code from Smyth: http://salty.oce.orst.edu/wave_analysis/SSF_index.html
% Smyth, W.D., J.N. Moum and J.D. Nash, 2011: ?Narrowband, high-frequency oscillations at the equator. Part II: Properties of shear instabilities", J. Phys. Oceanogr. 41, 412-428. 
% S.Tan, Scripps, 2019/04/15
%
% the example in this script is tested on tow-yo profiles, 
% only use DP.m and DP_0.m to process data
% use data from bottom to 200m above interface
% FG method, no background viscosity, diffusivity, 2nd scan
% S.Tan, Yantai, 2019/09/23

clear all
close all
clc

%% datapath & toolpath & etc. 
% main directory (change to your directory)
mdirec='/Users/tantanmeow/Desktop/WORK/2018-2019/Jesse/sp-tg/';

% data directory
datadirec='/Users/tantanmeow/Desktop/WORK/2018-2019/Jesse/proc_data/';
% add Matlab tools to path
addpath(genpath(strcat(mdirec,'tools/')));

%% parameters for the instability scan
tool = 2; % 1- FD ; 2- FG
% ----------  used in SSF.m or vTG.m  ----------  
% rangeing wave number for the first scan
% % 1) only scan over one direction (only one velocity component needed)
% k=-2:.01:-.5; 
% K=0; L=10.^(k); %successive values differing by a factor 10.^(.01)=1.0233
% disp(strcat('wave length ranges from ', num2str(2*pi/L(end)), 'm', 'to', num2str(2*pi/L(1)), 'm'))
% 2) only scan over two directions (both velocity component required)
k=[-1.5:.02:-.8]; K=[-fliplr(10.^(k)) 0 10.^(k)];
l=[-2:.02:-.8]; L=10.^(l); 
disp(strcat('x: wave length ranges from ', num2str(2*pi/K(end)), 'm', 'to', num2str(2*pi/10.^(k(1))), 'to inf'))
disp(strcat('y: wave length ranges from ', num2str(2*pi/L(end)), 'm', 'to', num2str(2*pi/L(1)), 'm'))
% x: wave length ranges from39.6442mto198.6918to inf
% y: wave length ranges from39.6442mto628.3185m
% (FD) background diffusion and dissipation
nu=0;
kap=0;
% (FD)
% boundary condition
% (1) velocity: 1=rigid (2nd order: w_{hat}=0; 4th order: w_{hat}=0, w_{hat}_z=0), 0=frictionless (2nd order: w_{hat}_z=kw_{hat}; 4th order: w_{hat}=0, w_{hat}_zz=0)
% (2) buoyancy: 1=insulating (b_{hat}_z=0 suggesting -\kappa\frac{\pb}{\pz}=0), 0=fixed-buoyancy (b_{hat}=0)                             
iBC1 = [1,0]; % at z=z(0)     
iBCN = iBC1;  % at z=z(N+1)

% (FG) or vertical and horizontal eddy viscosity, vertical and horizontal eddy diffusivity
Av=0;Ah=0;Kv=0;Kh=0;
% no need to specify BC : 
% Lian et al. (2020): This choice of basis functions is natural for impermeable, frictionless, constant-buoyancy
% boundaries, and may be adapted for insulating boundaries by expanding b_{hat} in terms of cosines
% rather than sines. But the method is not well-suited for rigid boundaries.

% ----------  used in DP.m  ----------  
% processed data grids (see HOWTO and BOT)
% (to note: a finer grid may better capture the critical level and you may want to use the new Matlab codes to do that since the old codes can be time consuming)
dz=10;
HOWTO = 1;% 1- linear interp; 2-segment mean
BOT=1; % data profile starts from 1-DEEPEST DATA POINT 2- REAL SEA BOTTOM

% ----------  used in FGM.m  ----------  
% bin-size of the measurments (resolution of data)
bs=8; % 10m or 16m? I'm not sure

% theredhold for k
k_thred = 1; % 1- 2.*pi./7/(bs/2); 2- do not apply

% perform 2nd scan
scan2 = 2;

%% below is an example for the sorted two-yo profiles
% load data
filename=dir(strcat(datadirec,'TY*'));
for filenum=1:length(filename)
    load(strcat(datadirec,filename(filenum).name))
    s = mlon; %number of profiles

II=nan(length(s),1);  % no. of zero-crossings
GR=nan(length(s),1);  % growth rate
CR=nan(length(s),1);  % real phase speed
CI=nan(length(s),1);  % imaginary phase speed
CL=nan(length(s),30); % critical level
KFGM=nan(length(s),2);% wave number (k, l) for the FGM 
botz=nan(length(s),1);% sea bottom depth (optional, data.hab required)
zw=[-6000:dz:0]';
W=nan(length(zw),length(s)); % fastest growing mode for w (complex)
% velocity profiles correspond to the FGM (wave vector parallel to flow)
V=nan(length(zw),length(s));Vz=nan(length(zw),length(s));Vzz=nan(length(zw),length(s));
% buoyancy, N squre and Richardson number profiles
B=nan(length(zw),length(s));N2=nan(length(zw),length(s));Ri=nan(length(zw),length(s));

GR_all=nan(length(s),length(K),length(L));  
% W_all=nan(length(zw),length(s),length(K),length(L)); 

vdata = v; udata = u; zdata = z; 
for profile=1:length(s)
    data.v=vdata(find(~isnan(vdata(:,profile))),profile);data.u=udata(find(~isnan(vdata(:,profile))),profile);data.z=zdata(find(~isnan(vdata(:,profile))),profile);
    data.hab=bdepth(profile)-max(depth(find(~isnan(vdata(:,profile))),profile));
    data.sigma4=sig4_sorted(find(~isnan(vdata(:,profile))),profile); %use sorted profiles that Jesse already generated
    data.b=b_sorted(find(~isnan(vdata(:,profile))),profile); 
    if isfield(data, 'hab')
        botz(profile)=data.z(end)-data.hab(end);
    end
    % D2=[];D1=[];% only keep data above bottom D1 meters to D2 meters
    D2=zo(profile)-botz(profile)+200;D1=0; %data from bottom to 200 m above the interface

    % scan the (k,l) plane once or twice to pick out FGM
    [v,vz,vzz,b,n2,Ri_r,zz,II(profile),GR(profile),CR(profile),CI(profile),WFGM,KFGM(profile,:),CL_FGM,II0,GR0,CR0,CI0,W0]=FGM_towyo_ex1(data,dz,D1,D2,HOWTO,BOT,K,L,nu,kap,Av,Ah,Kv,Kh,iBC1,iBCN,k_thred,bs,tool,scan2);
    CL(profile,1:length(CL_FGM))=CL_FGM;
    W(:,profile)=interp1(zz,WFGM,zw);V(:,profile)=interp1(zz,v,zw);Vz(:,profile)=interp1(zz,vz,zw);Vzz(:,profile)=interp1(zz,vzz,zw);
    B(:,profile)=interp1(zz,b,zw);N2(:,profile)=interp1(zz,n2,zw);Ri(:,profile)=interp1(zz,Ri_r,zw);

    % plot growth rates
    if ~isnan(GR(profile))
        figure(profile)
        SI_pic_1(GR0,K,L,bs,v,CR0,CI0,II0,GR(profile),CR(profile),CI(profile),KFGM(profile,:))
        print('-djpeg',[mdirec strcat('results/TG_SI_towyo/', filename(filenum).name(1:end-4), '_ex2_1_',num2str(profile))])
        % plot w
        figure(profile+50)
        [Ph]=SI_pic_2(WFGM,sqrt(KFGM(profile,1)^2+KFGM(profile,2)^2),zz,CL(profile,:));
        print('-djpeg',[mdirec strcat('results/TG_SI_towyo/', filename(filenum).name(1:end-4), '_ex2_2_',num2str(profile))])
    end
    % plot profiles
    figure(profile+100)
    DP_pic(v,vz,vzz,b,n2,Ri_r,zz)
    set(gcf, 'PaperPosition',[0 0 8 4]);
    print('-djpeg',[mdirec strcat('results/TG_SI_towyo/', filename(filenum).name(1:end-4), '_ex2_3_',num2str(profile))])

    close all
    
    GR_all(profile,:,:)=GR0;
%     W_all(:,profile,:,:)=interp1(zz,W0,zw);

end
LON = mlon; LAT = mlat;
eval(strcat('save TG_SI_',  filename(filenum).name(1:end-4), '_ex2.mat II GR CR CI CL KFGM K L zw W V Vz Vzz B N2 Ri LON LAT botz GR_all'))
end
