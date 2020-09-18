% The purpose of this analysis is to determine the phase speed of low mode
% waves using a 1D taylor goldstein analysis. 

% PARAMS
% pick and arbitrary profile to start with
pfl = 201;
idx = pfl;
% density of interface top
sig4i = 1045.93;  % double check this matches the choice in other analysis.
% Turbulent diffusivity (if constant)
Kv0 = 1e-2;
% To smooth profiles or not? Maybe diffusivity takes care of this?
dosmoothing = false;
% Smoothing amounts
bspan = 0.1;
vspan = 0.3;
% wavenumber for Taylor Goldstein
L = 100000;
k = 2*pi/L;
% Compute all modes (imode=1 gives fastest-growing unstable mode)
imode = 0;
% Skip some data to increase speed
step = 5;
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

% cut out some data for speed
u = u(1:step:end);
v = v(1:step:end);
sig4 = sig4(1:step:end);
b = b(1:step:end);
z = z(1:step:end);

% rotate the velocity into the depth mean direction
dz = 1*step;  % this is the bin size of each data point.
U = sum(u*dz)/sum(dz*use);
V = sum(v*dz)/sum(dz*use);
% need the - and 2*pi because atan2 returns values in the range [-pi, pi]
% and we want [0, 2*pi]
angle = -atan2(V, U) + 2*pi;
up = u*cos(angle) - v*sin(angle);
vp = u*sin(angle) + v*cos(angle);
% after doing the above, up should contain all the depth mean velocity.

%%%%%%% Taylor Goldstein analysis in the U direction

Mdiff = BaryL(z, 1, 6);  % This is the differentiation matrix.

% Smooth data
if dosmoothing
    us = smooth(z, up, vspan, 'rloess');
    bs = smooth(z, b, bspan, 'rloess');
    bz = Mdiff*bs;
else
    us = up;
    bs = b;
    bz = Mdiff*bs;
end

% Diffusivity profile
% kv = gamma*eps./bzs.^2;
[nz, ~] = size(bs);
Kv = Kv0*ones(nz, 1);
Kh = Kv;
Av = Kv;
Ah = Kv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Fourier-Galerkin method to compute growth rates & eigfns
% Compute Fourier integrals in advance
FG = vTG_FGprep(z, us, 0*us, bz, Av, Ah, Kv, Kh); 
% Compute growth rates & eigfns
[om, we, be] = vTG_FG(z, us, 0*us, Av, Ah, k, 0, imode, FG);
% read above as returning
% [frequency, w eigenvector, b eigenvector]

% phase speed
cp = -imag(om)/k;
% sort by phase speed
[cp, ind] = sort(cp,'ascend');
we = we(:, ind);
be = be(:, ind);
om = om(ind);
% horizontal velocity eigenfunction
ue = (1i/k)*Mdiff*we;

[~, imin] = min(cp);
[~, imax] = max(cp);
iplot = imin;

%%%%% ALL FIGURES BELOW!
% PARAMS
lw = 2;
% END PARAMS

figure
hold on
plot(lon, lat, 'k.')
plot(lon(idx), lat(idx), 'ro')

% figure
% plot(bs, z)
% 
% figure
% plot(bz, z)
% 
% figure
% plot(us, z)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the buoyancy and velocity profiles.
figure
subplot(2, 3, 1)
plot(bs, z, 'linewidth', lw)
title('b [m/s^2]')
ylabel('z [m]')
subplot(2, 3, 2)
hold on
plot(us, z, 'linewidth', lw)
title('u [m/s]')

% plot phase speeds
subplot(2, 3, 3)
hold on
plot(cp, 'bx')
xlabel('mode number')
ylabel('c_r [m/s]')
title('phase speeds')

hold on
plot(iplot, cp(iplot), 'rO', 'linewidth', 2, 'markersize',10)
axis tight;
yl = ylim;
ylim([yl(1)-.1, yl(2)+.1])

% Plot eigenfunctions for fastest phase speed.
subplot(2,3,4)
plot(real(we(:,iplot)),z,'linewidth',lw)
% hold on
% plot(imag(we(:,imode)),z,'r','linewidth',lw)
% legend('real','imag','location','southeast'); 
ylabel('z [m]')
xlabel('w eigfn')
title('fastest mode')

subplot(2,3,5)
plot(real(ue(:,iplot)), z, 'linewidth', lw)
% hold on
% plot(imag(ue(:,imode)),z,'r','linewidth',lw)
xlabel('u eigfn')
title(sprintf('c_r=%.3em/s', cp(iplot)))

subplot(2,3,6)
plot(real(be(:, iplot)), z, 'linewidth', lw)
% hold on
% plot(imag(be(:,imode)),z,'r','linewidth',lw)
xlabel('b eigfn')
title(sprintf('\\omega = (%.2e + %.2e i) s^{-1}', real(om(iplot)), imag(om(iplot))))
