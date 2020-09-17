%% Data processing step-2
% INPUTS:
% v,b,zz: (n,1) first element represents at the sea bottom or the deepest measurement
% dz is the grid size (usually = abs(mean(diff(zz))))
%
% OUTPUTS:
% v,vz,vzz,b_r,n2,Ri_r: fill upper layer flow reversals with stagnant and uniform fluid
%
% CALLS:
% ddz.m & ddz2.m
% turning_pt.m
%
% TO NOTICE:
% 
% S.Tan, IOCAS, 2020/09/17
function [v,vz,vzz,b,n2,Ri,zz]=DP_3(v,b,n2,zz,dz)

    l=find(isnan(v));
    v(l)=[];b(l)=[];n2(l)=[];zz(l)=[];
    
    Z=[zz(1):dz:0]';
    V=nan(size(Z));N2=nan(size(Z));B=nan(size(Z));
    V(1:length(zz))=v;N2(1:length(zz))=n2;
    
    % fill the upper-layer voids:
    l=find(isnan(N2));
    N2(l)=0;
    l=find(isnan(V));
    V(l)=0;
    % find the point of flow reversal
    [I,l]=turning_pt(V,0);
    % above that point, fill with stagnant and uniform fluid
    V(l(1)+1:end)=0;
    N2(l(1)+1:end)=0;
    
    B=cumtrapz(Z,N2)+b(1);% Integrate to get B

    % calculate vz, vzz, n2, Ri
    v = V; b = B; zz = Z;
    vz=ddz(zz)*v;
    vzz=ddz2(zz)*v;
    n2=ddz(zz)*b;
    Ri=n2./(vz).^2;
