%% This code plot pics on the results of shear instability
% INPUTS:
% w   : modes in respect of w
% K   : (ranging)
% i   : which K
% zz  : processed vertical grid
%
% OUTPUTS:
% pic of modes (w)
%
% TO NOTICE:
%
% S.Tan, Scripps, 2019/04/15
function [Ph]=SI_pic_2(w,K,zz,CL)
ms=15; lw=1.5;
subplot(1,2,1)
hold on
% h1=plot(real(w(:,1)),zz);
h1=plot(real(w(:,1)),zz,'linewidth',lw);
hold on
% h2=plot(imag(w(:,1)),zz,'-.');
h2=plot(imag(w(:,1)),zz,'-.','linewidth',lw);
hold on
plot(0,CL,'kp','markersize',ms,'linewidth',lw);
legend([h1,h2],'real','imag');
xlabel('Fatest growing mode w(z)');
ylabel('Depth');
title(strcat('k=',num2str(K)))
subplot(1,2,2)
hold on
% Ph=atan(imag(w(:,1))./real(w(:,1)));
% plot(Ph,(zz-mean(zz))./(2.*pi./K))%/pi*360 to degree %/180 dgree to rad
[Ph, Mag]=cart2pol(real(w(:,1)),imag(w(:,1)));
% plot(Ph,zz)
plot(Ph,zz,'linewidth',lw)
hold on
plot(0,CL,'kp','markersize',ms,'linewidth',lw);
xlabel('phase');
grid on
