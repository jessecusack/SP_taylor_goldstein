%% This is the main program to compute fastest growth rate of shear instabilities
%
% INPUTS:
% data: a Matlab structure that must contain fileds v, sigma4, z (optional: u, hab)
% v : the velocity
% sigma4 : potential density referenced to 4000 dbar
% hab : height above bottom in meters
% z : height in meters
% dz: z interval
% D1, D2: only keep data above bottom D1 meters to D2 meters (e.g., D2=1500;D1=0; )
% HOWTO: 1- linear interp; 2-segment mean
% BOT: data profile starts from 1-DEEPEST DATA POINT 2- REAL SEA BOTTOM
% K_init, L_init: initial guesses of wave numbers in x and y directions
% nu, kap: background viscosity, diffusivity
% iBC1, iBCN: boundary condition
% k_thred: theredhold for k ; 1- 2.*pi./7/(bs/2); 2- do not apply
% bs: bin-size of the measurments (resolution of data)
%
% OUTPUTS:
% KFGM & GR - wave number & fastest growth rate
%
% CALLS:
% FGM.m (change lines 45-50, 64-69, and 119-124 to process the profiles)
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
% the example in this script is tested on time-mean mporing profiles, 
% only use DP.m to process data
% assume rigid lid, no background viscosity, diffusivity
% S.Tan, IOCAS, 2019/09/06

clear all
close all
clc

%% datapath & toolpath & etc. 
% main directory (change to your directory)
mdirec='/Users/tantanmeow/Desktop/WORK/2018-2019/Jesse/sp-tg/';

% add Matlab tools to path
addpath(genpath(strcat(mdirec,'tools/')));

% load data
filename = 'mooring_averages.mat'; 
load(strcat(mdirec,'data/',filename))

%% parameters for the instability scan
% ----------  used in SSF.m  ----------  
% rangeing wave number for the first scan
% 1) only scan over one direction (only one velocity component needed)
k=-3:.01:0; 
K=0; L=10.^(k);
% % 2) only scan over two directions (both velocity component required)
% k=[-3:.1:0]; 
% K=[-fliplr(10.^(k)) 0 10.^(k)];L=10.^(k); 

% background diffusion and dissipation
nu=0;
kap=0;

% ----------  used in DP.m  ----------  
% processed data grids (see HOWTO and BOT)
% (to note: a finer grid may better capture the critical level and you may want to use the new Matlab codes to do that since the old codes can be time consuming)
dz=10;
HOWTO = 1;% 1- linear interp; 2-segment mean
BOT=1; % data profile starts from 1-DEEPEST DATA POINT 2- REAL SEA BOTTOM
D2=[];D1=[];% only keep data above bottom D1 meters to D2 meters

% boundary condition
% (1) velocity: 1=rigid, 0=frictionless
% (2) buoyancy: 1=insulating, 0=fixed-buoyancy
iBC1 = [1,0]; % at z=z(0)    
iBCN = iBC1;  % at z=z(N+1)

% ----------  used in FGM.m  ----------  
% bin-size of the measurments (resolution of data)
bs=8; % 10m or 16m? I'm not sure

% theredhold for k
k_thred = 1; % 1- 2.*pi./7/(bs/2); 2- do not apply

%% below is an example for the mean mooring profiles
s = {'P4','T8','T1','T2','T4','T6','T5', 'T12','M5','T9','T10','T11','P1','P3','T7'};

II=nan(length(s),1);  % no. of zero-crossings
GR=nan(length(s),1);  % growth rate
CR=nan(length(s),1);  % real phase speed
CI=nan(length(s),1);  % imaginary phase speed
CL=nan(length(s),10); % critical level
KFGM=nan(length(s),2);% wave number (k, l) for the FGM 
botz=nan(length(s),1);% sea bottom depth (optional, data.hab required)
zw=[-6000:dz:0]';
W=nan(length(zw),length(s)); % fastest growing mode for w (complex)
% velocity profiles correspond to the FGM (wave vector parallel to flow)
V=nan(length(zw),length(s));Vz=nan(length(zw),length(s));Vzz=nan(length(zw),length(s));
% buoyancy, N squre and Richardson number profiles
B=nan(length(zw),length(s));N2=nan(length(zw),length(s));Ri=nan(length(zw),length(s));

for profile=1:length(s)
    data=eval(s{profile});

    if isfield(data, 'hab')
        botz(profile)=data.z(end)-data.hab(end);
    end

    % scan the (k,l) plane and pick out FGM    
    [v,vz,vzz,b,n2,Ri_r,zz,II(profile),GR(profile),CR(profile),CI(profile),WFGM,KFGM(profile,:),CL_FGM,II0,GR0,CR0,CI0,W0]=FGM(data,dz,D1,D2,HOWTO,BOT,K,L,nu,kap,iBC1,iBCN,k_thred,bs);
    CL(profile,1:length(CL_FGM))=CL_FGM;
    W(:,profile)=interp1(zz,WFGM,zw);V(:,profile)=interp1(zz,v,zw);Vz(:,profile)=interp1(zz,vz,zw);Vzz(:,profile)=interp1(zz,vzz,zw);
    B(:,profile)=interp1(zz,b,zw);N2(:,profile)=interp1(zz,n2,zw);Ri(:,profile)=interp1(zz,Ri_r,zw);

    % plot growth rates
    if ~isnan(GR(profile))
        figure(profile)
        SI_pic_1(GR0,K,L,bs,v,CR0,CI0,II0,GR(profile),CR(profile),CI(profile),KFGM(profile,:))
        print('-djpeg',[mdirec strcat('results/TG_SI_mean/TG_SI_mean_ex1_1_',s{profile})])
        % plot w
        figure(profile+50)
        [Ph]=SI_pic_2(WFGM,sqrt(KFGM(profile,1)^2+KFGM(profile,2)^2),zz,CL(profile,:));
        print('-djpeg',[mdirec strcat('results/TG_SI_mean/TG_SI_mean_ex1_2_',s{profile})])
    end
    % plot profiles
    figure(profile+100)
    DP_pic(v,vz,vzz,b,n2,Ri_r,zz)
    set(gcf, 'PaperPosition',[0 0 8 4]);
    print('-djpeg',[mdirec strcat('results/TG_SI_mean/TG_SI_mean_ex1_3_',s{profile})])

    close all

end
site = s;
save TG_SI_mean_ex1.mat II GR CR CI CL KFGM K L zw W V Vz Vzz B N2 Ri site