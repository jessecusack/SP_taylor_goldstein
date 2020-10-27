function [u,p,K,SP,EF,cEF,BF,EFn,cEFn,eps,RH,LH,err]=FGM_energetics(z,U,k,l,nu,sig,w,b)

% INPUTS:
% z = vertical coordinate vector (evenly spaced)
% U = velocity profile
% k,l = wave vector 
% nu = viscosity
% sig = growth rate of FGM
% w = vertical velocity eigenfunction
% b = buoyancy eigenfunction
%
% OUTPUTS:
% u,p = horizontal velocity and pressure eigenfunction
% K = perturbation kinetic energy
% SP = shear production 
% EF = perturbation kinetic energy flux
% cEF = convergence of perturbation kinetic energy flux
% BF = buoyancy flux
% EFn = perturbation kinetic energy flux due to viscosity
% cEFn = (diffusion) convergence of perturbation kinetic energy flux due to viscosity
% eps =  (dissipation) viscous dissipation rate
% RH = right-hand side of kinetic energy budget; SP - EF_z + nu*Kzz - eps (idealy equal to sKE)
% LH = left-hand side of kinetic energy budget; evolution of the horizontally averaged perturbation kinetic energy: 2*real(sig)*K
% err = difference between RH and sKE, representing numerical errors
%
% CALLS:
% ddz.m & ddz2.m
% ref: Smyth and Carpenter (2019)
%
% TO NOTICE:
% assume vertical and horizontal viscosity are equal

ii = complex(0,1);
c = ii*sig/k;

% normalize base on w
for j=1:length(sig)
    cnorm=w(find(abs(w(:,j))==max(abs(w(:,j))),1),j);
    w(:,j)=w(:,j)/cnorm;
    b(:,j)=b(:,j)/cnorm;
end

% Compute pressure eigfn and horizontal velocity eigfn
% S.Tan 07/01/2019
Uz=ddz(z)*U;
Id=eye(length(z));
u=nan(size(w));p=nan(size(w));
for j=1:length(sig)
    u(:,j)=(sqrt(-1)/k).*ddz(z)*w(:,j);
    p(:,j)=-(sig(j)+sqrt(-1)*k*U).*(ddz(z)*w(:,j))./k^2 +(sqrt(-1)/k).*Uz.*w(:,j)+(nu/(k^2)).*((ddz2(z)-k^2.*Id)*(ddz(z)*w(:,j)));
    p(:,j)=(sqrt(-1)/k).*(Uz.*w(:,j)+(c(j)-U).*ddz(z)*w(:,j)-(sqrt(-1)/k*nu).*(ddz2(z)*ddz(z)*w(:,j)-k^2*ddz(z)*w(:,j)));
end

% Compute perturbation kinetic energy
wz=ddz(z)*w;
uz=ddz(z)*u;
K=(u.*conj(u)+w.*conj(w))/4;
% Compute shear production terms
uw=real(u.*conj(w))/2;
SP=-uw.*Uz;
% Compute perturbation kinetic energy flux
EF=real(p.*conj(w))/2;
cEF=-ddz(z)*EF;
% Compute buoyancy flux
BF=real(b.*conj(w))/2;
% Compute convergence of perturbation kinetic energy flux due to viscosity
EFn=-nu.*ddz(z)*K;
cEFn=nu.*ddz2(z)*K;
% Compute dissipation
eps=nu.*(uz.*conj(uz)+wz.*conj(wz)+(4*k^2).*K)./2;
for j=1:length(sig)
    LH(:,j)=2.*(real(sig(j)).*K(:,j));
end
RH=SP+BF+cEF+cEFn-eps;
err=abs(mean(RH-LH))./( mean(abs(RH).^2)+mean(abs(LH).^2) ).^.5;

return