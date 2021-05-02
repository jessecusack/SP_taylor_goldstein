%% Data processing step-0
% INPUTS:
% v,vz,vzz,b,n2,Ri,zz: (n,1) first element represents at the sea bottom or the deepest measurement
%
% OUTPUTS:
% v,vz,vzz,b_r,n2,Ri_r: get rid of N2<0
% 
% CALLS:
% ddz.m
%
% TO NOTICE:
% 
% S.Tan, IOCAS, 2020/09/17
% only extrapolate N2=nan instead of getting rid of N2<0
% S.Tan, IOCAS, 2021/01/19
function [v,vz,vzz,b,n2,Ri,zz]=DP_0(v,vz,vzz,b,n2,Ri,zz)
     l=find(isnan(v));
     v(l)=[];vz(l)=[];vzz(l)=[];b(l)=[];n2(l)=[];Ri(l)=[];zz(l)=[];
    
    % get rid of N2<0 !!!!!!!!! (this is a MUST, or GR increases with k, convective instability)
%     n2(find(n2<0))=nan;
    l=find(~isnan(n2));
    n2=interp1(zz(l),n2(l),zz,'pchip',999);
    l=find(n2==999);
    v(l)=[];vz(l)=[];vzz(l)=[];b(l)=[];n2(l)=[];Ri(l)=[];zz(l)=[];
        
%     b=cumtrapz(zz,n2)+b(1);% Integrate to get b
%     n2=ddz(zz)*b;
%     Ri=n2./(vz).^2;
%     l=find(n2<0);
%     if ismember(1,l)==1
%         n2(1)=mean(n2(2));
%     end
%     if ismember(length(n2),l)==1
%         n2(end)=mean(n2(end-1));
%     end
%     b=cumtrapz(zz,n2)+b(1);% Integrate to get b
%     n2=ddz(zz)*b;
%     Ri=n2./(vz).^2;
