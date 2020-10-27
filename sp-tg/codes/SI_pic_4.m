%% This code plot pics on the results of shear instability
% INPUTS:
%
% OUTPUTS:
% pic of buoyancy variances
%
% TO NOTICE:
%
% S.Tan, Yantai, 2020/10/02
function SI_pic_4(Bbud,kap,zz,CL_FGM)
fs=16;
lw=1.6;
ms=14;
    
subplot(161)
plot(Bbud.BV,zz,'k','linewidth',lw);
title('(a) BV','fontsize',fs,'fontweight','normal')
% plot(2*(imag(cph_h2{i}).*k)*Kbud.K,z_H,'k','linewidth',1.5);
% title('(a)2\sigma_r K','fontsize',12,'fontweight','normal')
ylabel('Depth(m)','fontsize',fs);
set(gca,'fontsize',fs);
grid on
text(min(Bbud.BV),zz(1)/5,strcat('\nu=\kappa='),'fontsize',fs);
text(min(Bbud.BV),zz(1)/4,strcat(num2str(kap)),'fontsize',fs);
subplot(162)
plot(Bbud.BVP,zz,'k','linewidth',lw);
title('(b) BVP','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
subplot(163)
plot(Bbud.BVF,zz,'k','linewidth',lw);
title('(c) BVF','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
subplot(164)
plot(Bbud.cBVF,zz,'k','linewidth',lw);
title('(d) BVF_z','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
subplot(165)
plot(Bbud.chi,zz,'k','linewidth',lw);
title('(e) \chi','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
subplot(166)
l(1)=plot(Bbud.LH_B,zz,'color',[.7 .7 .7],'linewidth',lw+2);
hold on
l(2)=plot(Bbud.RH_B,zz,'r','linewidth',lw);
title('(f) budget','fontsize',fs,'fontweight','normal')
set(gca,'fontsize',fs);
set(gca,'yticklabel',[]);
grid on
hold on
for i=1:length(CL_FGM)
    hold on
    plot([min([Bbud.LH_B;Bbud.RH_B]) max([Bbud.LH_B;Bbud.RH_B])],[CL_FGM(i) CL_FGM(i)],'k-.','linewidth',lw)
end
h=legend(l,'\sigma_r BV','BVP+BVF_z-\chi');

return