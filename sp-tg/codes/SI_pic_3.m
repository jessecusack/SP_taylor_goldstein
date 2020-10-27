%% This code plot pics on the results of shear instability
% INPUTS:
%
% OUTPUTS:
% pic of energetics
%
% TO NOTICE:
%
% S.Tan, Yantai, 2020/10/02
function SI_pic_3(Kbud,nu,zz,CL_FGM)
fs=16;
lw=1.6;
ms=14;

subplot(171)
plot(Kbud.K,zz,'k','linewidth',lw);
title('(a) K','fontsize',fs,'fontweight','normal')
% plot(2*(imag(cph_h2{i}).*k)*Kbud.K,z_H,'k','linewidth',1.5);
% title('(a)2\sigma_r K','fontsize',12,'fontweight','normal')
ylabel('Depth(m)','fontsize',fs)
set(gca,'fontsize',fs);
grid on
text(min(Kbud.K),zz(1)/5,strcat('\nu=\kappa='),'fontsize',fs);
text(min(Kbud.K),zz(1)/4,strcat(num2str(nu)),'fontsize',fs);
subplot(172)
plot(Kbud.SP,zz,'k','linewidth',lw);
title('(b) SP','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
subplot(173)
plot(Kbud.EF,zz,'k','linewidth',lw);
title('(c) EF','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
subplot(174)
plot(Kbud.cEF,zz,'k','linewidth',lw);
title('(d) EF_z','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
subplot(175)
plot(Kbud.BF,zz,'k','linewidth',lw);
title('(e) BF','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
subplot(176)
plot(Kbud.eps,zz,'k','linewidth',lw);
title('(f) \epsilon','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
subplot(177)
l(1)=plot(Kbud.LH,zz,'color',[.7 .7 .7],'linewidth',lw+2);
hold on
l(2)=plot(Kbud.RH,zz,'r','linewidth',lw);
title('(g) budget','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
hold on
for i=1:length(CL_FGM)
    hold on
    plot([min([Kbud.LH;Kbud.RH]) max([Kbud.LH;Kbud.RH])],[CL_FGM(i) CL_FGM(i)],'k-.','linewidth',1.5)
end
h=legend(l,'2\sigma_r K','SP-EF_z+BF-\epsilon');

return