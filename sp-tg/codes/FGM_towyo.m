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
% 4) second scan: reduce the grid spacing in the k, l plane; 
%    !!!!!!!!!!!  make sure to test with extended/reduced the domain,
%                 reduced the vertical grid spacing, FD or FG methods,
%                 inviscid vs viscous,
% S.Tan, IOCAS, 2019/09/22
% modified to suit TG_SI_towyo_ex4.m: 
% data process follows: DP, DP0, then low pass filter
% results include GR_all, CL_all (critical layer depth), Errs for perturbation Kinetic
% energy and perturbation density budgets
% S.Tan, Yantai, 2020/10/26
% add FD4 FD6 Chebyshev methods, delete codes for the second scan
% S.Tan, Yantai, 2020/10/28

function [v,vz,vzz,b,n2,Ri,zz,I_FGM,GR_FGM,CR_FGM,CI_FGM,W_FGM,K_FGM,CL_FGM,ERR_K_FGM,ERR_B_FGM,II0,GR0,CR0,CI0,W0,CL0,ERR_K0,ERR_B0]=FGM_towyo(data,dz,D1,D2,HOWTO,BOT,K_init,L_init,nu,kap,Av,Ah,Kv,Kh,iBC1,iBCN,k_thred,bs,tool,scan2,fl)

    tic
    K = K_init;L = L_init;  
    % data processing step-0
    [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz,D1,D2,K(1),L(1),HOWTO,BOT);
    [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz);
    % ------ first scan ------  
    KT0=nan(length(K),length(L));
    II0=nan(length(K),length(L));
    GR0=nan(length(K),length(L));
    CR0=nan(length(K),length(L));
    CI0=nan(length(K),length(L));
    W0=nan(length(zz),length(K),length(L));
    CL0=nan(length(K),length(L));
    ERR_K0=nan(length(K),length(L));
    ERR_B0=nan(length(K),length(L));
    for i=1:length(K)
        for j=1:length(L)           
            [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz,D1,D2,K(i),L(j),HOWTO,BOT);
            [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz); 
            % low-pass filter
            if fl>0
                v=butterworth_lp(v,fl/dz,4,length(v));
                vz=ddz(zz)*v;vzz=ddz2(zz)*v;
            end      
             % skip for min(Ri_r)>1/4 
             l=find(abs(Ri)<1/4);
             if ~isempty(l)
                 kt=sqrt(K(i)^2+L(j)^2);
                 if tool == 1
                    [sigs,w,bh]=SSF(zz,v,b,kt,0,nu,kap,iBC1,iBCN,1);
                    % perturbation kinetic energy
                    [u,p,Kbud.K,Kbud.SP,Kbud.EF,Kbud.cEF,Kbud.BF,Kbud.EFn,Kbud.cEFn,Kbud.eps,Kbud.RH,Kbud.LH,Kbud.err_K]=FGM_energetics(zz,v,kt,0,nu,sigs,w,bh);
                    % buoyancy variance
                    [Bbud.BV,Bbud.BVP,Bbud.LH_B,Bbud.RH_B,Bbud.chi,Bbud.BVF,Bbud.cBVF,Bbud.err_B]=FGM_buoyancy(zz,n2,kt,0,kap,sigs,w,bh);
                 elseif tool == 2
                     n2=BaryL(zz,1,6)*b; % differentiate buoyancy - 6th-order finite difference
                     vz=BaryL(zz,1,6)*v;
                     if length(zz)>7
                         vzz=BaryL(zz,2,6)*v;
                     else
                         vzz=nan(size(vz));
                     end
                     if length(Av)>1
                        FG = vTG_FGprep(zz,v,v*0,n2,Av,Ah,Kv,Kh); 
                        [sigs,w,bh]=vTG_FG(zz,v,v*0,Av,Ah,kt,0,1,FG);
                     else
                        FG = vTG_FGprep(zz,v,v*0,n2,ones(size(zz))*Av,ones(size(zz))*Ah,ones(size(zz))*Kv,ones(size(zz))*Kh); 
                        [sigs,w,bh]=vTG_FG(zz,v,v*0,ones(size(zz))*Av,ones(size(zz))*Ah,kt,0,1,FG);
                     end
                     % perturbation kinetic energy
                     [u,p,Kbud.K,Kbud.SP,Kbud.EF,Kbud.cEF,Kbud.BF,Kbud.EFn,Kbud.cEFn,Kbud.eps,Kbud.RH,Kbud.LH,Kbud.err_K]=FGM_energetics(zz,v,kt,0,Av,sigs,w,bh);
                     % buoyancy variance
                     [Bbud.BV,Bbud.BVP,Bbud.LH_B,Bbud.RH_B,Bbud.chi,Bbud.BVF,Bbud.cBVF,Bbud.err_B]=FGM_buoyancy(zz,n2,kt,0,Kv,sigs,w,bh);
                 elseif tool == 3 % new code FD-2
                     Div2.z = ddzv2(zz); Div2.zz = ddz2v2(zz); Div2.zzz = Div2.zz*Div2.z; Div2.zzzz = ddz4v2(zz);
                     n2=Div2.z*b; 
                     vz=Div2.z*v; 
                     vzz=Div2.zz*v; 
                     if length(Av)>1
                         [sigs,w,bh,e,u,p,Kbud,Bbud]=ParShear_AhAvKhKv_new_BS(zz,v,n2,Av,Ah,Kv,Kh,kt,0,Div2,iBC1,iBCN,1,1);
                     else
                         [sigs,w,bh,e,u,p,Kbud,Bbud]=ParShear_AhAvKhKv_new_BS(zz,v,n2,ones(size(zz))*Av,ones(size(zz))*Ah,ones(size(zz))*Kv,ones(size(zz))*Kh,kt,0,Div2,iBC1,iBCN,1,1);
                     end
                 elseif tool == 4 % FD-4
                      Div4.z = ddzv4(zz); Div4.zz = ddz2v4(zz); Div4.zzz = Div4.zz*Div4.z; Div4.zzzz = ddz4v4(zz);
                      n2=Div4.z*b; 
                      vz=Div4.z*v; 
                      vzz=Div4.zz*v; 
                      if length(Av)>1
                          [sigs,w,bh,e,u,p,Kbud,Bbud]=ParShear_AhAvKhKv_new_BS(zz,v,n2,Av,Ah,Kv,Kh,kt,0,Div4,iBC1,iBCN,1,1);
                      else
                          [sigs,w,bh,e,u,p,Kbud,Bbud]=ParShear_AhAvKhKv_new_BS(zz,v,n2,ones(size(zz))*Av,ones(size(zz))*Ah,ones(size(zz))*Kv,ones(size(zz))*Kh,kt,0,Div4,iBC1,iBCN,1,1);
                      end
                 elseif tool == 5 % FD-6
                      Div6.z = ddzv4(zz); Div6.zz = ddz2v4(zz); Div6.zzz = Div6.zz*Div6.z; Div6.zzzz = ddz4v6(zz);
                      n2=Div6.z*b; 
                      vz=Div6.z*v; 
                      vzz=Div6.zz*v; 
                      if length(Av)>1
                          [sigs,w,bh,e,u,p,Kbud,Bbud]=ParShear_AhAvKhKv_new_BS(zz,v,n2,Av,Ah,Kv,Kh,kt,0,Div6,iBC1,iBCN,1,1);
                      else
                          [sigs,w,bh,e,u,p,Kbud,Bbud]=ParShear_AhAvKhKv_new_BS(zz,v,n2,ones(size(zz))*Av,ones(size(zz))*Ah,ones(size(zz))*Kv,ones(size(zz))*Kh,kt,0,Div6,iBC1,iBCN,1,1);
                      end                   
                 elseif tool == 6 % Chebyshev
                 else
                    error('TOOL MUST BE 1 - 6')
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
                 ERR_K0(i,j)=Kbud.err_K;
                 ERR_B0(i,j)=Bbud.err_B;
                 CL0(i,j)=zz(find(Kbud.K==max(Kbud.K)));
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

    if ~isempty(l1)
        if scan2==1
            K = K_init;L = L_init;  
            II=II0;GR=GR0;CR=CR0;CI=CI0;W=W0;CL=CL0;ERR_K=ERR_K0;ERR_B=ERR_B0;
            [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz,D1,D2,K(l1),L(l2),HOWTO,BOT);
            [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz); 
            % low-pass filter
            if fl>0
                v=butterworth_lp(v,fl/dz,4,length(v));
            end
             if tool == 1
                 vz=ddz(zz)*v;vzz=ddz2(zz)*v;
             elseif tool == 2
                 n2=BaryL(zz,1,6)*b; % differentiate buoyancy - 6th-order finite difference
                 vz=BaryL(zz,1,6)*v;
                 if length(zz)>7
                     vzz=BaryL(zz,2,6)*v;
                 else
                     vzz=nan(size(vz));
                 end
             elseif tool == 3
                 Div2.z = ddzv2(zz); Div2.zz = ddz2v2(zz); Div2.zzz = Div2.zz*Div2.z; Div2.zzzz = ddz4v2(zz);
                 n2=Div2.z*b; 
                 vz=Div2.z*v; 
                 vzz=Div2.zz*v; 
             elseif tool == 4
                  Div4.z = ddzv4(zz); Div4.zz = ddz2v4(zz); Div4.zzz = Div4.zz*Div4.z; Div4.zzzz = ddz4v4(zz);
                  n2=Div4.z*b; 
                  vz=Div4.z*v; 
                  vzz=Div4.zz*v; 
             elseif tool == 5
                  Div6.z = ddzv4(zz); Div6.zz = ddz2v4(zz); Div6.zzz = Div6.zz*Div6.z; Div6.zzzz = ddz4v6(zz);
                  n2=Div6.z*b; 
                  vz=Div6.z*v; 
                  vzz=Div6.zz*v; 
             elseif tool == 6
                 error('fill in Chebyshev derivatives')
             else
                 error('TOOL MUST BE 1 - 6')
             end
        else
            error('only one scan available')
        end

        I_FGM = II(l1,l2);
        GR_FGM = GR(l1,l2);
        CR_FGM = CR(l1,l2);
        CI_FGM = CI(l1,l2);
        W_FGM = W(:,l1,l2);
        K_FGM = [K(l1) L(l2)];
        CL_FGM = CL(l1,l2);
        ERR_K_FGM = ERR_K(l1,l2);
        ERR_B_FGM = ERR_B(l1,l2);
    else
        I_FGM = nan; GR_FGM = nan; CR_FGM = nan; CI_FGM = nan; CL_FGM = nan; W_FGM = nan(size(v)); K_FGM = [nan nan];
        ERR_K_FGM = nan; ERR_B_FGM = nan;
        v = nan(size(v)); vz = nan(size(v)); vzz = nan(size(v));b = nan(size(v)); n2 = nan(size(v)); Ri = nan(size(v));
    end     
    toc
    
    
       