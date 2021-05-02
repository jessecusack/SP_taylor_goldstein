%% This code pre-processes data
% INPUTS:
% data: a Matlab structure that must contain fileds v, sigma4, z (optional: u, hab)
% dz: processed data grids
% D1, D2: only keep data above bottom D1 meters to D2 meters (e.g., D2=1500;D1=0; )
% K, L: wave number and direction of the wave vector
% HOWTO: 1- linear interp; 2-segment mean
% BOT: data profile starts from 1-DEEPEST DATA POINT 2- REAL SEA BOTTOM
%
% OUTPUTS:
% v,vz,vzz,b_r,n2,Ri_r: unaugmented profiles 
%
% CALLS:
% ddz.m & ddz2.m
%
% TO NOTICE:
% original codes
% S.Tan, Scripps, 2019/04/15
% speed rotated to the direction of wave
% S.Tan, IOCAS, 2020/09/06
% "data" better contains mean-state profiles 
% e.g., time averages, DP_1.m: sorted profiles 
% DP_2.m: artifically exclude overturns, BBLs 
% DP_3.m: fill upper layer flow reversals with stagnant and uniform fluid
% S.Tan, IOCAS, 2020/09/17
function [v,vz,vzz,b,n2,Ri,zz]=DP(data,dz,D1,D2,theta,HOWTO,BOT)
    
    v=data.v;sigma4=data.sigma4;z=data.z;b=nan(size(sigma4));
    if isfield(data, 'hab')
        hab=data.hab;
    else
        hab=nan(size(sigma4));
    end
    if isfield(data, 'u')
        u=data.u;
    else
        u=nan(size(sigma4));
    end
    if isfield(data, 'b')
        b=data.b;
    else
        % Convert density to buoyancy.
        g=9.81;
        l=find(~isnan(sigma4));
        b(l)=-g*(sigma4(l)-nanmean(sigma4))/nanmean(sigma4);
    end
    
    % z be negative
    if mean(z)>0
        z = -z;
    end
    % z be (n,1)
    if length(z(:,1))<2
        z = z';v = v';sigma4 = sigma4';b = b';hab = hab';
    end
    % z direct from bottom (z(1)) towards top (z(end))
    if mean(diff(z))<0
        z = flipud(z);v = flipud(v);sigma4 = flipud(sigma4);b = flipud(b);hab = flipud(hab);
    end

    % Evenly grid 
    l=find(~isnan(v));
    if BOT==2
        zz=[ceil(z(l(1))-hab(1)):dz:0]';
    elseif BOT==1
        zz=[ceil(z(l(1))):dz:0]';
    else
        error('BOT MUST BE 1 OR 2')   
    end      
        
    % Evenly grid for u, v, b
    V=nan(size(zz));U=nan(size(zz));B=nan(size(zz));
    if HOWTO==1
        l=find(~isnan(v));
        if length(l)>2
            V=interp1(z(l),v(l),zz,'linear');
        end
        l=find(~isnan(u));
        if length(l)>2
            U=interp1(z(l),u(l),zz,'linear');
        end
        l=find(~isnan(b));
        if length(l)>2  
            B=interp1(z(l),b(l),zz,'linear');
        end
    elseif HOWTO==2
        for i=2:length(zz)-1
            l=find(z>=zz(i)-dz/2&z<zz(i)+dz/2);
            V(i)=sum(v(l))./length(l);
            U(i)=sum(u(l))./length(l);
            B(i)=sum(b(l))./length(l);
        end
        l=find(z>=zz(end)-dz/2);
        V(end)=sum(v(l))./length(l);
        U(end)=sum(u(l))./length(l);
        B(end)=sum(b(l))./length(l);
        l=find(z<zz(1)+dz/2);
        V(1)=sum(v(l))./length(l);
        U(1)=sum(u(l))./length(l);
        B(1)=sum(b(l))./length(l);
    else 
        error('HOWTO MUST BE 1 OR 2')   
    end
    
    % near bottom (especially BOT==2): extrapolate V and use uniform B
%     l=find(isnan(V(1:round(length(find(~isnan(V)))/4))));ll=find(~isnan(V));
%     V(l)=interp1(zz(ll),V(ll),zz(l),'linear','extrap');
%     l=find(isnan(B(1:round(length(find(~isnan(B)))/4))));
%     B(l)=B(ll(1));
    % identify near bottom voids
    l=find(isnan(V));dl=find(diff(l)>1);
    if ~isempty(dl)
        ll=dl(1);l=l(1:ll);ll=find(~isnan(V));
    else
        l=find(isnan(V(1:round(length(find(~isnan(V)))/4))));ll=find(~isnan(V));
    end
    V(l)=interp1(zz(ll),V(ll),zz(l),'linear','extrap');
    l=find(isnan(B));dl=find(diff(l)>1);
    if ~isempty(dl)
        ll=dl(1);l=l(1:ll);ll=find(~isnan(B));
    else
        l=find(isnan(B(1:round(length(find(~isnan(B)))/4))));
    end
    B(l)=B(ll(1));    

    % Only keep data above seabed (last CTD point) D1 to D2
    v=V;u=U;b=B;
    if ~isempty(D2)
        ref=zz(1)+D2;
        za=zz-ref;
        l=find(abs(za)==min(abs(za)));
        v(l(1):end)=[];u(l(1):end)=[];b(l(1):end)=[];zz(l(1):end)=[];
    end
    if ~isempty(D1)
        ref=zz(1)+D1;
        za=zz-ref;
        l=find(abs(za)==min(abs(za)));
        v(1:l(1))=[];u(1:l(1))=[];b(1:l(1))=[];zz(1:l(1))=[];
    end

    % roated speed
    if isfield(data, 'u')
        spd = u.*cos(pi*theta/180) + v.*sin(pi*theta/180);
%         spd=(K.*u+L.*v)./sqrt(K.^2+L.^2);
    else
        spd=v;
    end

    % calculate unaugmented vz, vzz, n2, Ri_r for roated speed
    n2=nan(size(b));vz=nan(size(spd));vzz=nan(size(spd));
    l=find(~isnan(b));
    n2(l)=ddz(zz(l))*b(l);
    l=find(~isnan(spd));
    vz(l)=ddz(zz(l))*spd(l);
    vzz(l)=ddz2(zz(l))*spd(l);
    Ri=n2./(vz).^2;
    v=spd;
    
end

