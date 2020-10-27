%% Data processing step-2
% INPUTS:
% v,b,zz: (n,1) first element represents at the sea bottom or the deepest measurement
%
% OUTPUTS:
% v,vz,vzz,b_r,n2,Ri_r:  artifically exclude overturns, BBLs 
% 
% CALLS:
% ddz.m & ddz2.m
% butterworth_lp.m
%
% TO NOTICE:
% 1) low-pass filter for n2 and v is optional:
% fl is the band-width in meters
% dz is the grid size (usually = abs(mean(diff(zz))))
% 2) the thredhold (a small N2) for the BBL is arbitrarily, simply change in line 31
% S.Tan, IOCAS, 2020/09/17
function [v,vz,vzz,b,n2,Ri,zz]=DP_2(v,b,n2,zz,dz,fl)
    
    l=find(isnan(v));
    v(l)=[];b(l)=[];n2(l)=[];zz(l)=[];

    % If infos for low-pass filter are given
    % ---------- 1) filter n2 and v ---------- 
    if nargin<=5
        vz=ddz(zz)*v;vzz=ddz2(zz)*v;
    elseif nargin==6
        n2=butterworth_lp(n2,fl/dz,4,length(n2));
        v=butterworth_lp(v,fl/dz,4,length(v));
        vz=ddz(zz)*v;vzz=ddz2(zz)*v;
    end
    
    % ---------- 2) get rid of BBL & N2<0 ---------- 
    n2_bbl = 1e-6;
    % get rid of BBL
    l=find(n2<=n2_bbl);
    ll=diff(l);
    if ~isempty(ll)
    if max(ll)==1 %only upper or lower layer being well mixed
        l2=find(l<=length(zz)/4);
        n2(1:l(l2(end)))=[];b(1:l(l2(end)))=[];zz(1:l(l2(end)))=[];v(1:l(l2(end)))=[];vz(1:l(l2(end)))=[];vzz(1:l(l2(end)))=[];
    else
        ll1=find(ll>1);
        for i=1:length(ll1)
            if mean(n2(l(ll1(i))+1:l(ll1(i)+1)-1))<2*n2_bbl
                n2(l(ll1(i))+1:l(ll1(i)+1)-1)=n2_bbl;
            end
        end
        l=find(n2<=n2_bbl);
        ll=diff(l);ll1=find(ll>1);
        if ~isempty(ll1)
            n2(1:l(ll1(1))+1)=[];b(1:l(ll1(1))+1)=[];zz(1:l(ll1(1))+1)=[];v(1:l(ll1(1))+1)=[];vz(1:l(ll1(1))+1)=[];vzz(1:l(ll1(1))+1)=[];
        else
            n2(1:l(end)+1)=[];b(1:l(end)+1)=[];zz(1:l(end)+1)=[];v(1:l(end)+1)=[];vz(1:l(end)+1)=[];vzz(1:l(end)+1)=[];
        end
    end
    end
    
    % get rid of N2<0
    n2(find(n2<0))=nan;
    l=find(~isnan(n2));
    n2=interp1(zz(l),n2(l),zz,'pchip');
    % Integrate to get b
%     b=cumsum([n2(1)*dz/2;n2(2:end-1).*dz;n2(end)*dz/2])+b(1);
%     b=ddz(zz)\n2;
    b=cumtrapz(zz,n2)+b(1);
    n2=ddz(zz)*b;
    Ri=n2./(vz).^2;

