% Compute phase velocities and 1st mode eigenfunctions from Jon Nash's
% Columbia Plume data.
clear
close all
lw=2;
fs=16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the mean profiles of z, velocity and buoyancy gradient
% 1. Load the data
dat=load('Nash_data.txt');

% 2. Extract the data
V=dat(:,1);
rho=dat(:,2)*1000;
z=dat(:,3);

rho_0=mean(rho); 
% (You can choose the characteristic density value in different ways; 
% it doesn't make much difference to the results.)

% Convert density to buoyancy.
g=9.81;
B=-g*(rho-rho_0)/rho_0;
Bz=BaryL(z,1,6)*B; % differentiate buoyancy

% Plot the buoyancy and velocity profiles.
figure
subplot(2,3,1)
plot(B,z,'linewidth',2)
title('B [m/s^2]')
ylabel('z [m]')
subplot(2,3,2)
plot(V,z,'linewidth',2)
title('V [m/s]')

% Specify viscosity, diffusivity
nu=1.e-6; % molecular value for water
Av=nu*ones(size(z));Ah=Av;
Kv=Av/7;Kh=Kv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Fourier-Galerkin method to compute growth rates & eigfns
% First specify wavenumber. (For the hydrostatic limit, just choose a small, 
% but nonzero, value).
k=.00001;
% Compute all modes (imode=1 gives fastest-growing unstable mode)
imode=0;
% Compute Fourier integrals in advance
FG = vTG_FGprep(z,V,V*0,Bz,Av,Ah,Kv,Kh); 
% Compute growth rates & eigfns
[sigs,w,b]=vTG_FG(z,V,V*0,Av,Ah,k,0,imode,FG);

% Convert growth rates to phase speeds. 
% (For instability problems, work with the growth rate. 
% For wave problems, the phase velocity is more useful.)
cphs=-imag(sigs)/k;
% Sort modes by phase speed, fastest first.
[cphs,ind]=sort(cphs,'descend');
w=w(:,ind);b=b(:,ind);
sigs=sigs(ind);
% horizontal velocity eigenfunction
u=(sqrt(-1)/k)*BaryL(z,1,6)*w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
% plot phase speeds
subplot(2,3,3)
plot(cphs,'bx')
xlabel('mode number')
ylabel('c_r [m/s]')
title('phase speeds')
imode=find(abs(cphs)==max(abs(cphs)),1);
hold on
plot(imode,cphs(imode),'rO','linewidth',2,'markersize',10)
axis tight;
yl=ylim;
ylim([yl(1)-.1,yl(2)+.1])

% Plot eigenfunctions for fastest phase speed.
subplot(2,3,4)
plot(real(w(:,imode)),z,'linewidth',lw)
hold on
plot(imag(w(:,imode)),z,'r','linewidth',lw)
legend('real','imag','location','southeast'); 
ylabel('z [m]')
xlabel('w eigfn')
title('fastest mode')

subplot(2,3,5)
plot(real(u(:,imode)),z,'linewidth',lw)
hold on
plot(imag(u(:,imode)),z,'r','linewidth',lw)
xlabel('u eigfn')
title(sprintf('c_r=%.3em/s',cphs(imode)))

subplot(2,3,6)
plot(real(b(:,imode)),z,'b','linewidth',lw)
hold on
plot(imag(b(:,imode)),z,'r','linewidth',lw)
xlabel('b eigfn')
title(sprintf('\\sigma=(%.2e,%.2e)s^{-1}',real(sigs(imode)),imag(sigs(imode))))