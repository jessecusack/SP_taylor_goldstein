%% Data processing step-1
% INPUTS:
% v,b,zz: (n,1) first element represents at the sea bottom or the deepest measurement
%
% OUTPUTS:
% v,vz,vzz,b_r,n2,Ri_r: sorted profiles (buoyant waters on top)
% 
% CALLS:
% ddz.m & ddz2.m
%
% TO NOTICE:
% 
% S.Tan, IOCAS, 2020/09/17
function [v,vz,vzz,b,n2,Ri,zz]=DP_1(v,b,zz)

    l=find(isnan(v));
    v(l)=[];b(l)=[];zz(l)=[];
    
    % only sort the buoyancy profile
    [bs,order]=sort(b);
    b=bs;

    % calculate vz, vzz, n2, Ri
    n2=nan(size(b));vz=nan(size(v));vzz=nan(size(v));
    l=find(~isnan(b));
    n2(l)=ddz(zz(l))*b(l);
    l=find(~isnan(v));
    vz(l)=ddz(zz(l))*v(l);
    vzz(l)=ddz2(zz(l))*v(l);
    Ri=n2./(vz).^2;
