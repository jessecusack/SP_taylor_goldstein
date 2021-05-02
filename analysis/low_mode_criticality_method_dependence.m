% The purpose of this analysis is to determine the phase speed of low mode
% waves using a 1D taylor goldstein analysis. 

% PARAMS
% pick and arbitrary profile to start with
pfl = 212;
idx = pfl;
% density of interface top
sig4i = 1045.93; %1045.93;  % double check this matches the choice in other analysis.
% Turbulent diffusivity (if constant)
Kv0 = 1e-2;
% To smooth profiles or not? Maybe diffusivity takes care of this?
dosmoothing = true;
% Smoothing amounts
zlp = 50;  % low pass wavelength in m
zd = 1;  % sampling wavelength in m
Np = zlp/zd;  % number of data points per wavelength
mlp = 1/zlp;  % filter wavenumber cpm
ms = 1/zd;  % sampling wavenumber cpm
% wavenumber for Taylor Goldstein
L = [1000 10000 100000 1000000];
% Compute all modes (imode=1 gives fastest-growing unstable mode)
imode = 0;
% Skip some data to increase speed
step = 5;
% Plotting params
lw = 2;
iplot = 1;
% Boundary conditions
iBC1 = [0 0];
iBCN = [0 0];
% END PARAMS

% LOAD DATA
% This data file is generated by the script 'stack_towyos.py'.
data_file = "../proc_data/stacked_towyos.nc";
u_ = ncread(data_file, "u");
v_ = ncread(data_file, "v");
b_ = ncread(data_file, "b_sorted");
z_ = ncread(data_file, "z");
sig4_ = ncread(data_file, "sig4_sorted");
lon = ncread(data_file, "lon");
lat = ncread(data_file, "lat");

% clean up data and remove data above the interface
u = u_(idx, :);
v = v_(idx, :);
sig4 = sig4_(idx, :);
b = b_(idx, :);

use = ~isnan(b) & ~isnan(u) & ~isnan(v) & (sig4 > sig4i);

% also make everything a column vector... erugh matlab
u = u(use)';
v = v(use)';
sig4 = sig4(use)';
b = b(use)';
z = z_(use);

U = sum(u*zd)/sum(zd*use);
V = sum(v*zd)/sum(zd*use);
% need the - and 2*pi because atan2 returns values in the range [-pi, pi]
% and we want [0, 2*pi]
angle = -atan2(V, U) + 2*pi;
up = u*cos(angle) - v*sin(angle);
vp = u*sin(angle) + v*cos(angle);
% after doing the above, up should contain all the depth mean velocity.

% Diffusivity profile
% kv = gamma*eps./bzs.^2;
Kv = Kv0*ones(size(b));

%% %%%%% Taylor Goldstein analysis in the U direction

Mdiff = BaryL(z, 1, 6);  % This is the differentiation matrix.

% Smooth data
if dosmoothing
    disp('Doing lowpass')
    us = lowpass(up, mlp, ms, 'ImpulseResponse','iir');
    bs = lowpass(b, mlp, ms, 'ImpulseResponse','iir');
    sig4s = lowpass(sig4, mlp, ms, 'ImpulseResponse','iir');
    bz = Mdiff*bs;
    zs = z;
    Kvs = Kv;
else
    us = up;
    bs = b;
    sig4s = sig4;
    bz = Mdiff*bs;
    zs = z;
    Kvs = Kv;
    
end


% cut out some data for speed
us = us(1:step:end);
% v = v(1:step:end);
sig4s = sig4s(1:step:end);
bs = bs(1:step:end);
bz = bz(1:step:end);
zs = zs(1:step:end);
Kvs = Kvs(1:step:end);

% Other diffusivity/viscosity params
Khs = Kvs;
Avs = Kvs;
Ahs = Kvs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Fourier-Galerkin method to compute growth rates & eigfns
% Compute Fourier integrals in advance
FG = vTG_FGprep(zs, us, 0*us, bz, Avs, Ahs, Kvs, Khs); 

figure
subplot(2, 3, 1);
hold on
plot(b, z, 'linewidth', lw)
plot(bs, zs, 'linewidth', lw)
title('b [m/s^2]')
ylabel('z [m]')

subplot(2, 3, 2);
hold on
plot(up, z, 'linewidth', lw)
plot(us, zs, 'linewidth', lw)
title('u [m/s]')

s3 = subplot(2, 3, 3);
hold on
xlabel('mode number')
ylabel('c_r [m/s]')
title('phase speeds')

s4 = subplot(2,3,4);
hold on
ylabel('z [m]')
xlabel('w eigfn')
title('fastest mode')

s5 = subplot(2,3,5);
hold on
xlabel('u eigfn')

for L_ = L
    k = 2*pi/L_;
    
    [om, we, be] = vTG_FG(zs, us, 0*us, Avs, Ahs, k, 0, imode, FG);
    
    % phase speed
    cp = -imag(om)/k;
    % sort by phase speed
    [cp, ind] = sort(cp,'ascend');
    we = we(:, ind);
    be = be(:, ind);
    om = om(ind);
    % horizontal velocity eigenfunction
    ue = (1i/k)*BaryL(zs, 1, 3)*we;
    
    plot(s3, cp, 'x')
    plot(s4, real(we(:, iplot))/max(abs(real(we(:, iplot)))), zs, 'linewidth', lw)
    plot(s5, real(ue(:, iplot))/max(abs(real(ue(:, iplot)))), zs, 'linewidth', lw)
    
end

legend(s3, '1000', '10000', '100000', '1000000')


%%%%%%%%%%%%%%%%%%%%%%
% Older method with finite differences

figure
subplot(2, 3, 1);
hold on
plot(b, z, 'linewidth', lw)
plot(bs, zs, 'linewidth', lw)
title('b [m/s^2]')
ylabel('z [m]')

subplot(2, 3, 2);
hold on
plot(up, z, 'linewidth', lw)
plot(us, zs, 'linewidth', lw)
title('u [m/s]')

s3 = subplot(2, 3, 3);
hold on
xlabel('mode number')
ylabel('c_r [m/s]')
title('phase speeds')

s4 = subplot(2,3,4);
hold on
ylabel('z [m]')
xlabel('w eigfn')
title('fastest mode')

s5 = subplot(2,3,5);
hold on
xlabel('u eigfn')

for L_ = L
    k = 2*pi/L_;
    
    [om, we, be] = SSF(zs, us, bs, k, 0, Kv0, Kv0, iBC1, iBCN, imode);
    
    % phase speed
    cp = -imag(om)/k;
    % sort by phase speed
    [cp, ind] = sort(cp,'ascend');
    we = we(:, ind);
    be = be(:, ind);
    om = om(ind);
    % horizontal velocity eigenfunction
    ue = (1i/k)*BaryL(zs, 1, 3)*we;
    
    plot(s3, cp, 'x')
    plot(s4, real(we(:, iplot))/max(abs(real(we(:, iplot)))), zs, 'linewidth', lw)
    plot(s5, real(ue(:, iplot))/max(abs(real(ue(:, iplot)))), zs, 'linewidth', lw)
    
end

legend(s3, '1000', '10000', '100000', '1000000')