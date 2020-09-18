%% This code plot pics on the results of shear instability
% INPUTS:
% GR  : growth rate (rad/hr)
% K   : (ranging)
% bs  : (bin size of data)
% v   : v raw
% CR  : real phase speed 
% CI  : imag phase speed
%
% OUTPUTS:
% pic of k - fastest growth rate
% pic of Howard's semicircle theorem
%
% TO NOTICE:
% sort according to 1) growth rate 2) gravest stable modes by the turing point
%
% S.Tan, Scripps, 2019/04/15

function SI_pic_1(GR,K,L,bs,v,CR,CI,II,GR_max,CR_max,CI_max,K_max)
COLOR=colormap(jet);
COLOR=COLOR(1:3:end,:);
fs = 14; lw = 2; ms = 15;
if length(K)>1&&length(L)>1
    subplot(211)
    contourf(K,L,GR')
    hold on
    plot([2.*pi./7/(bs/2) 2.*pi./7/(bs/2)],[L(1) L(end)],'K');
    hold on
    plot([-2.*pi./7/(bs/2) -2.*pi./7/(bs/2)],[L(1) L(end)],'K');
    hold on
    plot([K(1) K(end)],[2.*pi./7/(bs/2) 2.*pi./7/(bs/2)],'K');
    xlabel('wave number (k)')
    ylabel('wave number (l)')
    title('Growth rate');
    colormap(COLOR);
    bar=colorbar;
    hold on
    plot(K_max(1),K_max(2),'kp','markersize',ms,'linewidth',lw);
    text(K(1),L(end)/5-L(1),strcat('Fastest growing wave'),'fontsize',fs,'color','k');
    text(K(1),L(end)/3-L(1),strcat('Wave length=',num2str(2*pi/sqrt(K_max(1)^2+K_max(2)^2)),' m'),'fontsize',fs,'color','k');
    text(K(1),2*L(end)/3-L(1),strcat('Growth rate= 2\pi/',num2str(2*pi/GR_max/3600),'hours'),'fontsize',fs,'color','k');
    grid on
else
    subplot(211)
    if length(K)>1
        for i=1:length(K)
            hold on
            if II(i)>9
                plot(log10(K(i)),GR(i),'*','color',COLOR(end,:));
            else
                plot(log10(K(i)),GR(i),'*','color',COLOR(II(i)+1,:));
            end
        end
        text_idx = K(1);
    else
        for i=1:length(L)
            hold on
            if II(i)>9
                plot(log10(L(i)),GR(i),'*','color',COLOR(end,:));
            else
                plot(log10(L(i)),GR(i),'*','color',COLOR(II(i)+1,:));
            end
        end 
        text_idx = L(1);
    end
    hold on
    plot([log10(2.*pi./7/(bs/2)) log10(2.*pi./7/(bs/2))],[0 GR_max],'K');
    ylabel('growth rate')
    xlabel('log(wave number (k))')
    colormap(COLOR);
    bar=colorbar;
    caxis([0,11]);
    set(bar,'ylim',[0,11],'yticklabel',{'','0','1','2','3','4','5','6','7','8','9','>9'})
    hold on
    K_max = max(K_max);
    plot(log10(K_max),GR_max,'rp','markersize',ms,'linewidth',lw);
    text(log10(text_idx),GR_max,strcat('Fastest growing wave'),'fontsize',fs,'color','k');
    text(log10(text_idx),GR_max-2*GR_max/10,strcat('Wave length=',num2str(2*pi/K_max),' m'),'fontsize',fs,'color','k');
    text(log10(text_idx),GR_max-3*GR_max/10,strcat('Growth rate= 2\pi/',num2str(2*pi/GR_max/3600),'hours'),'fontsize',fs,'color','k');
    grid on  
end

% pic of Howard's semicircle theorem
subplot(212)
plot(CR,CI,'k*')
hold on
plot(CR_max,CI_max,'rp','markersize',ms,'linewidth',lw);
hold on
th = 0:pi/50:pi;
Umax=max(v);
Umin=min(v);
r=(Umax-Umin)/2;
xunit = r * cos(th) + (Umax+Umin)/2;
yunit = r * sin(th) + 0;
plot(xunit, yunit,'k');
ylabel('c_i')
xlabel('c_r')
set(gcf, 'color', 'w')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 6*96 3*96]);
set(gcf, 'InvertHardcopy', 'off')
grid on
hold on

