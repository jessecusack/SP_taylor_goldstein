%% This code plot processed profiles
% INPUTS:
% v,vz,vzz,b,n2,Ri,zz: (1,n) columns
%
% OUTPUTS:
%
% TO NOTICE:
%
% S.Tan, Scripps, 2019/04/09
function DP_pic(v,vz,vzz,b,n2,Ri,zz)
    ms=10; lw=1.5;
    subplot(161)
    hold on
    plot(v,zz,'linewidth',lw)
    title('V');
    subplot(162)
    hold on
    plot(vz,zz,'linewidth',lw)
    title('Vz');
    subplot(163)
    hold on
    plot(vzz,zz,'linewidth',lw)
    title('Vzz');
    subplot(164)
    hold on
    plot(b,zz,'linewidth',lw)
    title('B');
    subplot(165)
    hold on
    plot(n2,zz,'linewidth',lw)
    hold on
    plot(vz.^2./4,zz,'r-.','linewidth',lw)
    title('N2');
    subplot(166)
    hold on
    plot(log10(Ri),zz,'linewidth',lw)
    hold on
    Ref = log10(1/4).*ones(size(Ri));
    Ref(isnan(Ri))=nan;
    plot(Ref, zz,'k','linewidth',lw)
    l=find(Ri<0);
    hold on
    plot(zeros(size(zz(l))),zz(l),'kp','markersize',ms,'linewidth',lw)
    l=find(Ri<1/4);
    hold on
    plot(zeros(size(zz(l))),zz(l),'r^','markersize',ms,'linewidth',lw)
    title('Ri');    
