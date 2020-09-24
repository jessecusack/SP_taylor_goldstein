%% This code picks out the fastest growing mode 
% INPUTS:
% data: a Matlab structure that must contain fileds v, sigma4, z (optional: u, hab)
% dz: z interval
% D1, D2: only keep data above bottom D1 meters to D2 meters (e.g., D2=1500;D1=0; )
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
% 0) SSF.m was designed to handle latitudial background flows (U), 
%    lines 111-115 may be problematic if one use meridional flow speed (V) instead of U while using true k and l.
%    A easy solution is to use the flow speed that is projected onto the wave vector (our case v) 
%    and to use the wave-vector magnitude kt=sqrt(k^2+l^2) for k and 0 for l
% 1) wavelenght criteria: 
% Consider a hyperbolic tangent shear layer: U=tanh(z/h),
% the ratio of wavelength to shear layer thickness (h, the shear layer boundaries are roughly z=+-h) is about 7 (valid in the unstratified and stratified cases)
% 2) call DP_0.m: artificially get rid of N2<0, or GR increases with k - convective instability
% 3) skip scan when min(Ri_r)>1/4  
% S.Tan, IOCAS, 2019/09/06
% 4) haven't wrote the codes for the second scan: 
%   Inspriation from Smyth et al. (2011): 
%       i) halve the grid spacing in the k, l plane; 
%       ii) extend the domain; (call DP_3.m?) 
%       iii) reduce the vertical grid spacing to 1 m. (new Matlab tool is more efficient)
% S.Tan, IOCAS, 2019/09/17

function [v,vz,vzz,b,n2,Ri,zz,I_FGM,GR_FGM,CR_FGM,CI_FGM,W_FGM,K_FGM,CL_FGM,II,GR,CR,CI,W]=FGM(data,dz,D1,D2,HOWTO,BOT,K_init,L_init,nu,kap,Av,Ah,Kv,Kh,iBC1,iBCN,k_thred,bs,tool,scan2)

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
                 if tool == 1
                    [sigs,w]=SSF(zz,v,b,kt,0,nu,kap,iBC1,iBCN,1);
                 elseif tool == 2
                     n2=BaryL(zz,1,6)*b; % differentiate buoyancy - 6th-order finite difference
                     vz=BaryL(zz,1,6)*v;vzz=BaryL(zz,2,6)*v;
                     if length(Av)>1
                        FG = vTG_FGprep(zz,v,v*0,n2,Av,Ah,Kv,Kh); 
                        [sigs,w]=vTG_FG(zz,v,v*0,Av,Ah,kt,0,1,FG);
                     else
                        FG = vTG_FGprep(zz,v,v*0,n2,ones(size(zz))*Av,ones(size(zz))*Ah,ones(size(zz))*Kv,ones(size(zz))*Kh); 
                        [sigs,w]=vTG_FG(zz,v,v*0,Av,Ah,kt,0,1,FG);
                     end
                 else
                     error('TOOL MUST BE 1 OR 2')
                 end
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

    % seek for FGM
    % wavelength therdhold (Hazel, 1972?)
    if k_thred == 2
       [l1, l2]=find(GR0==max(max(GR0)));
    elseif k_thred == 1
        a=GR0; a(find(KT0>2*pi/7/(bs/2)))=0;
        [l1, l2]=find(GR0==max(max(a)));       
    else
       error('k_thred must be 1 or 2');
    end
    II=II0;GR=GR0;CR=CR0;CI=CI0;W=W0;

    if ~isempty(l1)
    if scan2==2
    % ------- second scan ------%  
    % 1. reduced grid spacing (1/4) in the k, l plane
        if length(K)==1
            K_iter = K;
        else
            K_iter = 10.^(log10(K(l1))-log10(K(2)/K(1))*3:log10(K(2)/K(1))/4:log10(K(l1))+log10(K(2)/K(1))*3);
        end
        if length(L)==1
            L_iter = L;
        else
            L_iter = 10.^(log10(L(l2))-log10(L(2)/L(1))*3:log10(L(2)/L(1))/4:log10(L(l2))+log10(L(2)/L(1))*3);
        end            
    % 2. reduce the vertical grid spacing to dz/2
        % data processing step-0
        [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz/2,D1,D2,K_iter(1),L_iter(1),HOWTO,BOT);
        % data processing step-1 (optional)
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_1(v,b,zz); 
        % data processing step-2
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_2(v,b,n2,zz,dz/2,100); % low pass filter :(... ,dz/2,fl)
        % data processing step-3
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_3(v,b,n2,zz,dz/2);
        [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz);
        KT1=nan(length(K_iter),length(L_iter));
        II1=nan(length(K_iter),length(L_iter));
        GR1=nan(length(K_iter),length(L_iter));
        CR1=nan(length(K_iter),length(L_iter));
        CI1=nan(length(K_iter),length(L_iter));
        W1=nan(length(zz),length(K_iter),length(L_iter));

        for i=1:length(K_iter)
            for j=1:length(L_iter)           
                [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz/2,D1,D2,K_iter(i),L_iter(j),HOWTO,BOT);
                % data processing step-1 (optional)
            %     [v,vz,vzz,b,n2,Ri,zz]=DP_1(v,b,zz); 
                % data processing step-2
            %     [v,vz,vzz,b,n2,Ri,zz]=DP_2(v,b,n2,zz,dz/2,100); % low pass filter :(... ,dz/2,fl)
                % data processing step-3
            %     [v,vz,vzz,b,n2,Ri,zz]=DP_3(v,b,n2,zz,dz/2);
                [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz); 

                 % skip for min(Ri_r)>1/4 
                 l=find(abs(Ri)<1/4);
                 if ~isempty(l)
                     kt=sqrt(K_iter(i)^2+L_iter(j)^2);
                     if tool == 1
                        [sigs,w]=SSF(zz,v,b,kt,0,nu,kap,iBC1,iBCN,1);
                     elseif tool == 2
                         n2=BaryL(zz,1,6)*b; % differentiate buoyancy - 6th-order finite difference
                         vz=BaryL(zz,1,6)*v;vzz=BaryL(zz,2,6)*v;
                         if length(Av)>1
                            FG = vTG_FGprep(zz,v,v*0,n2,Av,Ah,Kv,Kh); 
                            [sigs,w]=vTG_FG(zz,v,v*0,Av,Ah,kt,0,1,FG);
                         else
                            FG = vTG_FGprep(zz,v,v*0,n2,ones(size(zz))*Av,ones(size(zz))*Ah,ones(size(zz))*Kv,ones(size(zz))*Kh); 
                            [sigs,w]=vTG_FG(zz,v,v*0,Av,Ah,kt,0,1,FG);
                         end
                     else
                         error('TOOL MUST BE 1 OR 2')
                     end
                     inst=real(sigs)/kt;
                     cphs=-imag(sigs)/kt;
                     ww=real(w(:,:));
                     I=turning_pt(ww,0); % turing point          
                     II1(i,j)=I(1);
                     GR1(i,j)=real(sigs(1));
                     CR1(i,j)=cphs(1);
                     CI1(i,j)=inst(1);
                     W1(:,i,j)=w;
                     KT1(i,j)=kt;
                 end
            end
        end
    % 3. estimate the location of the GR maximum using a biquadratic fit 
        if k_thred == 2
           [l1, l2]=find(GR1==max(max(GR1)));
        elseif k_thred == 1
            a=GR1; a(find(KT1>2*pi/7/(bs/2)))=0;
            [l1, l2]=find(GR1==max(max(a)));       
        else
           error('k_thred must be 1 or 2');
        end
        K = K_iter;L = L_iter;  
        II=II1;GR=GR1;CR=CR1;CI=CI1;W=W1;
        [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz/2,D1,D2,K(l1),L(l2),HOWTO,BOT);
        % data processing step-1 (optional)
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_1(v,b,zz); 
        % data processing step-2
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_2(v,b,n2,zz,dz/2,100); % low pass filter :(... ,dz/2,fl)
        % data processing step-3
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_3(v,b,n2,zz,dz/2);
        [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz);   
         if tool == 2
             n2=BaryL(zz,1,6)*b; % differentiate buoyancy - 6th-order finite difference
             vz=BaryL(zz,1,6)*v;vzz=BaryL(zz,2,6)*v;
         end
    elseif scan2==1
        K = K_init;L = L_init;  
        II=II0;GR=GR0;CR=CR0;CI=CI0;W=W0; 
        [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz,D1,D2,K(l1),L(l2),HOWTO,BOT);
        % data processing step-1 (optional)
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_1(v,b,zz); 
        % data processing step-2
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_2(v,b,n2,zz,dz,100); % low pass filter :(... ,dz,fl)
        % data processing step-3
    %     [v,vz,vzz,b,n2,Ri,zz]=DP_3(v,b,n2,zz,dz);
        [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz);   
         if tool == 2
             n2=BaryL(zz,1,6)*b; % differentiate buoyancy - 6th-order finite difference
             vz=BaryL(zz,1,6)*v;vzz=BaryL(zz,2,6)*v;
         end
    else
        error('SCAN2 MUST BE 1 OR 2')
    end
    
        I_FGM = II(l1,l2);
        GR_FGM = GR(l1,l2);
        CR_FGM = CR(l1,l2);
        CI_FGM = CI(l1,l2);
        W_FGM = W(:,l1,l2);
        K_FGM = [K(l1) L(l2)];

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
    
    
       