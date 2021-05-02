%% This code picks out mode families and their fastest growing mode
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
% GR_lm: ignore growth rates that are shorter than this (hr^-1)
% imode: how many modes to retain
% 
% OUTPUTS:
% FGM properties of 10 modes families: GR_family,CR_family,CI_family,CL_family,ERR_family,K_family,L_family,Knetm_family,SPm_family,BFm_family,EFnm_family,cEFnm_family,epsm_family,RHm_family,LHm_family,
% ends_family: 1-FGM at the end of the wavenumer range; 0-not
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
% 1) first scan: for each kt and psi, find the mode families and pick out
% the FGM (retain imode numbers of the growing modes)
% 2) ignore growth rates that are shorter than GR_lm
% 3) plot the histogram of the critical levels of retained modes and pick 
% out peaks with numbers larger than 3
% 4) at each critical level family, separate the families based on psi, mode 
% families that contain more than three modes are chosen, growth rate of the 
% FGM are computed


% not support choosing FD method yet
% S.Tan, IOCAS, 2021/01/31
function [GR_family,CR_family,CI_family,CL_family,ERR_family,K_family,L_family,Knetm_family,SPm_family,BFm_family,EFnm_family,cEFnm_family,epsm_family,RHm_family,LHm_family,Knet_family,SP_family,BF_family,EFn_family,cEFn_family,eps_family,RH_family,LH_family,ends_family]=FGM_modefamily(data,dz,D1,D2,HOWTO,BOT,wavenumber,theta,Av,Ah,Kv,Kh,GR_lm,imode)

    tic
    K = wavenumber;L = theta;  
    % ------ 1) for each kt and psi, find the mode families and pick out the FGM ------  
    GR = nan(length(K),length(L),imode);
    CR = nan(length(K),length(L),imode);
    CI = nan(length(K),length(L),imode);
    CL = nan(length(K),length(L),imode);
    ERR = nan(length(K),length(L),imode);

    Knet = nan(length(K),length(L),imode);
    SP = nan(length(K),length(L),imode);
    BF = nan(length(K),length(L),imode);
    EFn = nan(length(K),length(L),imode);
    cEFn = nan(length(K),length(L),imode);
    eps = nan(length(K),length(L),imode);
    RH = nan(length(K),length(L),imode);
    LH = nan(length(K),length(L),imode);
    Knetm = nan(length(K),length(L),imode);
    SPm = nan(length(K),length(L),imode);
    BFm = nan(length(K),length(L),imode);
    EFnm = nan(length(K),length(L),imode);
    cEFnm = nan(length(K),length(L),imode);
    epsm = nan(length(K),length(L),imode);
    RHm = nan(length(K),length(L),imode);
    LHm = nan(length(K),length(L),imode);

    % compute critical levels for all the growing modes
    for i=1:length(K)
        for j=1:length(L)    
            [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz,D1,D2,L(j),HOWTO,BOT);
            [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz); 
            l=find(abs(Ri)<1/4|min(n2)<0);
            kt=K(i);

            FG = vTG_FGprep(zz,v,v*0,n2,Av,Ah,Kv,Kh); 
            [sigs,w,bh]=vTG_FG(zz,v,v*0,Av,Ah,kt,0,0,FG);

            index_gr = find(real(sigs).*3600>0);
            if length(index_gr)>4
                inst=real(sigs(index_gr))/kt;
                cphs=-imag(sigs(index_gr))/kt;
                gr = real(sigs(index_gr));
                cl = nan(size(cphs));
                err = nan(size(cphs));

                K0 = nan(size(cphs));
                SP0 = nan(size(cphs));
                BF0 = nan(size(cphs));
                EFn0 = nan(size(cphs));
                cEFn0 = nan(size(cphs));
                eps0 = nan(size(cphs));
                RH0 = nan(size(cphs));
                LH0 = nan(size(cphs));
                K0m = nan(size(cphs));
                SP0m = nan(size(cphs));
                BF0m = nan(size(cphs));
                EFn0m = nan(size(cphs));
                cEFn0m = nan(size(cphs));
                eps0m = nan(size(cphs));
                RH0m = nan(size(cphs));
                LH0m = nan(size(cphs));

                for n=1:length(index_gr)
                    [u,p,Kbud.K,Kbud.SP,Kbud.EF,Kbud.cEF,Kbud.BF,Kbud.EFn,Kbud.cEFn,Kbud.eps,Kbud.RH,Kbud.LH,Kbud.err_K]=FGM_energetics(zz,v,kt,Av,sigs(index_gr(n)),w(:,index_gr(n)),bh(:,index_gr(n)));
                    cl(n)=zz(find(Kbud.K==max(Kbud.K)));
                    err(n)=Kbud.err_K;
                    % at K = max(K)
                    K0(n)=Kbud.K(find(Kbud.K==max(Kbud.K)));
                    SP0(n)=Kbud.SP(find(Kbud.K==max(Kbud.K)));
                    BF0(n)=Kbud.BF(find(Kbud.K==max(Kbud.K)));
                    EFn0(n)=Kbud.EFn(find(Kbud.K==max(Kbud.K)));
                    cEFn0(n)=Kbud.cEFn(find(Kbud.K==max(Kbud.K)));
                    eps0(n)=Kbud.eps(find(Kbud.K==max(Kbud.K)));
                    RH0(n)=Kbud.RH(find(Kbud.K==max(Kbud.K)));
                    LH0(n)=Kbud.LH(find(Kbud.K==max(Kbud.K)));
                    % average within +-1std of K masked by max(K)/10
                    K_ = Kbud.K;
                    K_(find(abs(K_)<max(Kbud.K)/10)) = nan;
                    K0m(n)=mean(Kbud.K(find(Kbud.K>=nanstd(K_))));
                    SP0m(n)=mean(Kbud.SP(find(Kbud.K>=nanstd(K_))));
                    BF0m(n)=mean(Kbud.BF(find(Kbud.K>=nanstd(K_))));
                    EFn0m(n)=mean(Kbud.EFn(find(Kbud.K>=nanstd(K_))));
                    cEFn0m(n)=mean(Kbud.cEFn(find(Kbud.K>=nanstd(K_))));
                    eps0m(n)=mean(Kbud.eps(find(Kbud.K>=nanstd(K_))));
                    RH0m(n)=mean(Kbud.RH(find(Kbud.K>=nanstd(K_))));
                    LH0m(n)=mean(Kbud.LH(find(Kbud.K>=nanstd(K_))));
                end
                h=histogram(cl,round(length(cl)/2), 'Orientation', 'horizontal');
                bins = h.BinEdges;
                counts = h.Values;
                close

                % now find peaks
                [pks,locs] = findpeaks(counts);
                % indexes of the fastest growing mode in each mode family
                index_family = nan(size(pks));
                gr_family = nan(size(pks));
                for n=1:length(pks)
                    bins_ = bins(locs(n):locs(n)+1);
                    index_family_ = find(cl>=bins_(1) & cl<bins_(2));%probably wrong
                    if length(index_family_)<pks(n)
                        index_family_ = find(cl>=bins_(1) & cl<=bins_(2));
                    end
                    gr_family_ = gr(index_family_);

                    [gr_family(n),index_family(n)] = max(gr_family_);
                    index_family(n) = index_family_(index_family(n));
                end
                cr_family = cphs(index_family);
                ci_family = inst(index_family);
                cl_family = cl(index_family);
                err_family = err(index_family);
                w_family = w(:,index_family);

                K0_family = K0(index_family);
                SP0_family = SP0(index_family);
                BF0_family = BF0(index_family);
                EFn0_family = EFn0(index_family);
                cEFn0_family = cEFn0(index_family);
                eps0_family = eps0(index_family);
                RH0_family = RH0(index_family);
                LH0_family = LH0(index_family);
                K0m_family = K0m(index_family);
                SP0m_family = SP0m(index_family);
                BF0m_family = BF0m(index_family);
                EFn0m_family = EFn0m(index_family);
                cEFn0m_family = cEFn0m(index_family);
                eps0m_family = eps0m(index_family);
                RH0m_family = RH0m(index_family);
                LH0m_family = LH0m(index_family);

                [gr_family, index_family] = sort(gr_family, 'descend');
                cr_family = cr_family(index_family);
                ci_family = ci_family(index_family);
                cl_family = cl_family(index_family);
                err_family = err_family(index_family);
                w_family = w_family(:,index_family);

                K0_family = K0_family(index_family);
                SP0_family = SP0_family(index_family);
                BF0_family = BF0_family(index_family);
                EFn0_family = EFn0_family(index_family);
                cEFn0_family = cEFn0_family(index_family);
                eps0_family = eps0_family(index_family);
                RH0_family = RH0_family(index_family);
                LH0_family = LH0_family(index_family);
                K0m_family = K0m_family(index_family);
                SP0m_family = SP0m_family(index_family);
                BF0m_family = BF0m_family(index_family);
                EFn0m_family = EFn0m_family(index_family);
                cEFn0m_family = cEFn0m_family(index_family);
                eps0m_family = eps0m_family(index_family);
                RH0m_family = RH0m_family(index_family);
                LH0m_family = LH0m_family(index_family);

                if length(index_family)<=imode
                    GR(i,j,1:length(index_family)) = gr_family;
                    CR(i,j,1:length(index_family)) = cr_family;
                    CI(i,j,1:length(index_family)) = ci_family;
                    CL(i,j,1:length(index_family)) = cl_family;
                    ERR(i,j,1:length(index_family)) = err_family;
                    Knet(i,j,1:length(index_family)) = K0_family;
                    SP(i,j,1:length(index_family)) = SP0_family;
                    BF(i,j,1:length(index_family)) = BF0_family;
                    EFn(i,j,1:length(index_family)) = EFn0_family;
                    cEFn(i,j,1:length(index_family)) = cEFn0_family;
                    eps(i,j,1:length(index_family)) = eps0_family;
                    RH(i,j,1:length(index_family)) = RH0_family;
                    LH(i,j,1:length(index_family)) = LH0_family;
                    Knetm(i,j,1:length(index_family)) = K0m_family;
                    SPm(i,j,1:length(index_family)) = SP0m_family;
                    BFm(i,j,1:length(index_family)) = BF0m_family;
                    EFnm(i,j,1:length(index_family)) = EFn0m_family;
                    cEFnm(i,j,1:length(index_family)) = cEFn0m_family;
                    epsm(i,j,1:length(index_family)) = eps0m_family;
                    RHm(i,j,1:length(index_family)) = RH0m_family;
                    LHm(i,j,1:length(index_family)) = LH0m_family;
                else
                    GR(i,j,:) = gr_family(1:imode);
                    CR(i,j,:) = cr_family(1:imode);
                    CI(i,j,:) = ci_family(1:imode);
                    CL(i,j,:) = cl_family(1:imode);
                    ERR(i,j,:) = err_family(1:imode);
                    Knet(i,j,:) = K0_family(1:imode);
                    SP(i,j,:) = SP0_family(1:imode);
                    BF(i,j,:) = BF0_family(1:imode);
                    EFn(i,j,:) = EFn0_family(1:imode);
                    cEFn(i,j,:) = cEFn0_family(1:imode);
                    eps(i,j,:) = eps0_family(1:imode);
                    RH(i,j,:) = RH0_family(1:imode);
                    LH(i,j,:) = LH0_family(1:imode);
                    Knetm(i,j,:) = K0m_family(1:imode);
                    SPm(i,j,:) = SP0m_family(1:imode);
                    BFm(i,j,:) = BF0m_family(1:imode);
                    EFnm(i,j,:) = EFn0m_family(1:imode);
                    cEFnm(i,j,:) = cEFn0m_family(1:imode);
                    epsm(i,j,:) = eps0m_family(1:imode);
                    RHm(i,j,:) = RH0m_family(1:imode);
                    LHm(i,j,:) = LH0m_family(1:imode);
                end
            end
        end
    end

    % ------ 2) ignore growth rates that are shorter than GR_lm ------ 
    index_all = find(GR.*3600>=GR_lm);
    K_all=nan(size(GR));
    L_all=nan(size(GR));
    for j=1:length(L)
        for n=1:imode
            K_all(:,j,n)=K;
        end
    end
    for i=1:length(K)
        for n=1:imode
            L_all(i,:,n)=L;
        end
    end
    GR_all = GR(index_all);
    CR_all = CR(index_all);
    CI_all = CI(index_all);
    CL_all = CL(index_all);
    ERR_all = ERR(index_all);
    K_all = K_all(index_all);
    L_all = L_all(index_all);

    Knet_all = Knet(index_all);
    SP_all = SP(index_all);
    BF_all = BF(index_all);
    EFn_all = EFn(index_all);
    cEFn_all = cEFn(index_all);
    eps_all = eps(index_all);
    RH_all = RH(index_all);
    LH_all = LH(index_all);

    Knetm_all = Knetm(index_all);
    SPm_all = SPm(index_all);
    BFm_all = BFm(index_all);
    EFnm_all = EFnm(index_all);
    cEFnm_all = cEFnm(index_all);
    epsm_all = epsm(index_all);
    RHm_all = RHm(index_all);
    LHm_all = LHm(index_all);

    % ------ 3) histogram of the critical level retained modes and pick out peaks ------ 
    h=histogram(CL_all, zz(end)-zz(1), 'Orientation', 'horizontal');
    bins = h.BinEdges;
    counts = h.Values;
    close
    % one mode family contains more than three modes
    [pks,locs] = findpeaks(counts,'MinPeakProminence',3);

    % ------ 4) separate the families based on psi and pick out the FGM ------ 
    % pick out mode families from different psi
    GR_family = nan(length(pks),10);
    CR_family = nan(length(pks),10);
    CI_family = nan(length(pks),10);
    CL_family = nan(length(pks),10);
    ERR_family = nan(length(pks),10);
    K_family = nan(length(pks),10);
    L_family = nan(length(pks),10);

    Knet_family = nan(length(pks),10);
    SP_family = nan(length(pks),10);
    BF_family = nan(length(pks),10);
    EFn_family = nan(length(pks),10);
    cEFn_family = nan(length(pks),10);
    eps_family = nan(length(pks),10);
    RH_family = nan(length(pks),10);
    LH_family = nan(length(pks),10);
    Knetm_family = nan(length(pks),10);
    SPm_family = nan(length(pks),10);
    BFm_family = nan(length(pks),10);
    EFnm_family = nan(length(pks),10);
    cEFnm_family = nan(length(pks),10);
    epsm_family = nan(length(pks),10);
    RHm_family = nan(length(pks),10);
    LHm_family = nan(length(pks),10);
    ends_family = zeros(length(pks),10);

    for n=1:length(pks)
        bins_ = bins(locs(n):locs(n)+1);
        index_family_ = find(CL_all>=bins_(1) & CL_all<bins_(2));%probably wrong
        if length(index_family_)<pks(n)
            index_family_ = find(CL_all>=bins_(1) & CL_all<=bins_(2));
        end
        GR_family_ = GR_all(index_family_);
        CR_family_ = CR_all(index_family_);
        CI_family_ = CI_all(index_family_);
        CL_family_ = CL_all(index_family_);
        ERR_family_ = ERR_all(index_family_);
        K_family_ = K_all(index_family_);
        L_family_ = L_all(index_family_);

        Knet_family_ = Knet_all(index_family_);
        SP_family_ = SP_all(index_family_);
        BF_family_ = BF_all(index_family_);
        EFn_family_ = EFn_all(index_family_);
        cEFn_family_ = cEFn_all(index_family_);
        eps_family_ = eps_all(index_family_);
        RH_family_ = RH_all(index_family_);
        LH_family_ = LH_all(index_family_);

        Knetm_family_ = Knetm_all(index_family_);
        SPm_family_ = SPm_all(index_family_);
        BFm_family_ = BFm_all(index_family_);
        EFnm_family_ = EFnm_all(index_family_);
        cEFnm_family_ = cEFnm_all(index_family_);
        epsm_family_ = epsm_all(index_family_);
        RHm_family_ = RHm_all(index_family_);
        LHm_family_ = LHm_all(index_family_);

    %     figure
    %     subplot(221)
    %     plot(K_family_,GR_family_.*3600,'.')
    %     subplot(223)
    %     plot(K_family_,CL_family_,'.')
    %     subplot(224)
    %     plot(L_family_,CL_family_,'.')
    %     subplot(222)
    %     plot(K_family_,L_family_,'.')

%         figure
        count = 0;
        for j=1:length(L)
            ll = find(L_family_==L(j));
            if ~isempty(ll)
                % 2) one mode family contains more than three modes
                if length(ll)>=3
                    count = count + 1;
                    [K_family_pick,pick] = sort(K_family_(ll));       
                    index_family = ll(pick);
                    [GR_family(n,count),index_max] = max(GR_family_(index_family));
                    CR_family(n,count) = CR_family_(index_family(index_max));
                    CI_family(n,count) = CI_family_(index_family(index_max));
                    CL_family(n,count) = CL_family_(index_family(index_max));
                    ERR_family(n,count) = ERR_family_(index_family(index_max));
                    K_family(n,count) = K_family_(index_family(index_max));
                    L_family(n,count) = L_family_(index_family(index_max));
                    
                    Knet_family(n,count) = Knet_family_(index_family(index_max));
                    SP_family(n,count) = SP_family_(index_family(index_max));
                    BF_family(n,count) = BF_family_(index_family(index_max));
                    EFn_family(n,count) = EFn_family_(index_family(index_max));
                    cEFn_family(n,count) = cEFn_family_(index_family(index_max));
                    eps_family(n,count) = eps_family_(index_family(index_max));
                    RH_family(n,count) = RH_family_(index_family(index_max));
                    LH_family(n,count) = LH_family_(index_family(index_max));

                    Knetm_family(n,count) = Knetm_family_(index_family(index_max));
                    SPm_family(n,count) = SPm_family_(index_family(index_max));
                    BFm_family(n,count) = BFm_family_(index_family(index_max));
                    EFnm_family(n,count) = EFnm_family_(index_family(index_max));
                    cEFnm_family(n,count) = cEFnm_family_(index_family(index_max));
                    epsm_family(n,count) = epsm_family_(index_family(index_max));
                    RHm_family(n,count) = RHm_family_(index_family(index_max));
                    LHm_family(n,count) = LHm_family_(index_family(index_max));
                    
                    % if the fastest growing mode is at the ends of the wave length range
                    if index_max==1||index_max==length(index_family)
                        ends_family(n,count) = 1;   
                    else
                        ends_family(n,count) = 0;   
                    end

%                     subplot(221)
%                     hold on
%                     plot(K_family_(index_family),GR_family_(index_family).*3600,'o-')
%                     hold on
%                     plot(K_family_(index_family(index_max)),GR_family_(index_family(index_max)).*3600,'rp','linewidth',1.5);
%                     xlabel('k');ylabel('gr')    
%                     subplot(222)
%                     hold on
%                     plot(K_family_(index_family),L_family_(index_family),'o-')
%                     hold on
%                     plot(K_family_(index_family(index_max)),L_family_(index_family(index_max)),'rp','linewidth',1.5);
%                     xlabel('k');ylabel('psi')    
%                     subplot(223)
%                     hold on
%                     plot(K_family_(index_family),CL_family_(index_family),'o-')
%                     hold on
%                     plot(K_family_(index_family(index_max)),CL_family_(index_family(index_max)),'rp','linewidth',1.5);
%                     xlabel('k');ylabel('cl')    
%                     subplot(224)
%                     hold on
%                     plot(L_family_(index_family),CL_family_(index_family),'o-')
%                     hold on
%                     plot(L_family_(index_family(index_max)),CL_family_(index_family(index_max)),'rp','linewidth',1.5);
%                     xlabel('psi');ylabel('cl')                
                end
            end
        end
        [a,ll_] = nanmax(GR_family(n,:));
        if ~isempty(ll_)
            ll = find(L_family_==L_family(n,ll_)); 
            [K_family_pick,pick] = sort(K_family_(ll));       
            index_family = ll(pick);
            hold on
            plot(2.*pi./K_family_(index_family),GR_family_(index_family).*3600,'o-')
            hold on
            plot(2.*pi./K_family(n,ll_),GR_family(n,ll_).*3600,'rp','linewidth',1.5);
            xlabel('wave length (m)');ylabel('growth rate (/hr)')    
        end
    end
    toc
