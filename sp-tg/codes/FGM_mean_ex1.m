%% This code picks out the fastest growing mode 
% INPUTS:
% data: a Matlab structure that must contain fileds v, sigma4, z (optional: u, hab)
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
% v,vz,vzz,b,n2,Ri,zz: profiles used to compute the FGM
% I_FGM,GR_FGM,CR_FGM,CI_FGM,W_FGM,K_FGM,CL_FGM: FGM
% 
% CALLS:
% SSF.m
% DP.m: basic data processing
% DP_0.m: get rid of N2<0
% optional:
% DP_1.m: sort buoyancy profile
% DP_2.m: low-pass filter, artifically exclude overturns, BBLs 
% DP_3.m: fill upper layer flow reversals with stagnant and uniform fluid
%
% TO NOTICE:
% modified FGM.m for TG_SI_mean_ex1.m
% S.Tan, IOCAS, 2019/09/18

function [v,vz,vzz,b,n2,Ri,zz,I_FGM,GR_FGM,CR_FGM,CI_FGM,W_FGM,K_FGM,CL_FGM,II0,GR0,CR0,CI0,W0]=FGM_mean_ex1(data,dz,D1,D2,HOWTO,BOT,K_init,L_init,nu,kap,iBC1,iBCN,k_thred,bs)

    tic
    K = K_init;L = L_init;  
    % data processing step-0
    [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz,D1,D2,K(1),L(1),HOWTO,BOT);
    % data processing step-1 (optional)
%     [v,vz,vzz,b,n2,Ri,zz]=DP_1(v,b,zz); 
    % data processing step-2
%     [v,vz,vzz,b,n2,Ri,zz]=DP_2(v,b,n2,zz,dz,100); % low pass filter :(... ,dz,fl)
    % data processing step-3
%     [v,vz,vzz,b,n2,Ri,zz]=DP_3(v,b,n2,zz,dz);
    [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz);

    % ------ first scan ------  
    KT0=nan(length(K),length(L));
    II0=nan(length(K),length(L));
    GR0=nan(length(K),length(L));
    CR0=nan(length(K),length(L));
    CI0=nan(length(K),length(L));
    W0=nan(length(zz),length(K),length(L));
    for i=1:length(K)
        for j=1:length(L)           
            [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz,D1,D2,K(i),L(j),HOWTO,BOT);
            % data processing step-1 (optional)
        %     [v,vz,vzz,b,n2,Ri,zz]=DP_1(v,b,zz); 
            % data processing step-2
        %     [v,vz,vzz,b,n2,Ri,zz]=DP_2(v,b,n2,zz,dz,100); % low pass filter :(... ,dz,fl)
            % data processing step-3
        %     [v,vz,vzz,b,n2,Ri,zz]=DP_3(v,b,n2,zz,dz);
            [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz); 
            
             % skip for min(Ri_r)>1/4 
             l=find(abs(Ri)<1/4);
             if ~isempty(l)
                 kt=sqrt(K(i)^2+L(j)^2);
                 [sigs,w]=SSF(zz,v,b,kt,0,nu,kap,iBC1,iBCN,1);
                 inst=real(sigs)/kt;
                 cphs=-imag(sigs)/kt;
                 ww=real(w(:,:));
                 I=turning_pt(ww,0); % turing point          
                 II0(i,j)=I(1);
                 GR0(i,j)=real(sigs(1));
                 CR0(i,j)=cphs(1);
                 CI0(i,j)=inst(1);
                 W0(:,i,j)=w;
                 KT0(i,j)=kt;
             end
        end
    end

    % ------ seek for FGM ------ 
    % wavelength therdhold (Hazel, 1972?)
    if k_thred == 2
       [l1, l2]=find(GR0==max(max(GR0)));
    elseif k_thred == 1
        a=GR0; a(find(KT0>2*pi/7/(bs/2)))=0;
        [l1, l2]=find(GR0==max(max(a)));       
    else
       error('k_thred must be 1 or 2');
    end
    
    % ------- second scan ------%  
    %                           %
    %       to be added soon    %
    %                           %
    % --------------------------%  


    II=II0;GR=GR0;CR=CR0;CI=CI0;W=W0;
    if ~isempty(l1)
        I_FGM = II(l1,l2);
        GR_FGM = GR(l1,l2);
        CR_FGM = CR(l1,l2);
        CI_FGM = CI(l1,l2);
        W_FGM = W(:,l1,l2);
        K_FGM = [K(l1) L(l2)];
        [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz,D1,D2,K(l1),L(l2),HOWTO,BOT);
        % data processing step-1 (optional)
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_1(v,b,zz); 
        % data processing step-2
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_2(v,b,n2,zz,dz,100); % low pass filter :(... ,dz,fl)
        % data processing step-3
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_3(v,b,n2,zz,dz);
        [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz);   
        
        % critical level
        [tp,l]=turning_pt(v,CR_FGM);CL_FGM=[];
        if tp == 0
           CL_FGM = nan;
        else
           for i=1:tp
               CL_FGM(i) = interp1([v(l(i)) v(l(i)+1)],[zz(l(i)) zz(l(i)+1)],CR_FGM);
           end
        end
    else
        I_FGM = nan; GR_FGM = nan; CR_FGM = nan; CI_FGM = nan; CL_FGM = nan; W_FGM = nan(size(v)); K_FGM = [nan nan];
        v = nan(size(v)); vz = nan(size(v)); vzz = nan(size(v));b = nan(size(v)); n2 = nan(size(v)); Ri = nan(size(v));
    end 
    
    toc
    
    
       