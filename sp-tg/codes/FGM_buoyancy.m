%% this code compute thhe buoyancy budget
function [BV,BVP,LH_B,RH_B,chi,BVF,cBVF,err_B]=FGM_buoyancy(z,Bz,k,l,kap,sig,w,b)

% INPUTS:
% z = vertical coordinate vector (evenly spaced)
% Bz = buoyancy gradient
% k,l = wave vector 
% kap = diffusivity
% sig = growth rate of FGM
% w = vertical velocity eigenfunction
% b = buoyancy eigenfunction
%
% OUTPUTS:
% BV = buoyancy variance
% BVP = buoyancy variance production;
% LH_B = left-hand side of buoyancy variance budget;
% RH_B = right-hand side of buoyancy variance budget;
% chi = buoyancy variance dissipation rate
% BVF = flux
% cBVF = flux convergence;
% err_B = relative error;
%
% CALLS:
% ddz.m & ddz2.m
% ref: Smyth and Carpenter (2019)
%
% TO NOTICE:
% assume vertical and horizontal diffusivity are equal

ii = complex(0,1);
c = ii*sig/k;

% normalize base on w
for j=1:length(sig)
    cnorm=w(find(abs(w(:,j))==max(abs(w(:,j))),1),j);
    w(:,j)=w(:,j)/cnorm;
    b(:,j)=b(:,j)/cnorm;
end

% Buoyancy variance budget
bz=ddz(z)*b;

BV=(b.*conj(b))/2;
LH_B=( diag(real(sig))*BV' )';

BVP=-diag(Bz)*real(conj(b).*w)/2;
chi=kap*abs(bz).^2/2+kap*k*abs(b).^2/2;

BVF=-kap*real(bz.*conj(b))/2;
cBVF=ddz(z)*BVF;
RH_B=BVP+cBVF-chi;
err_B= mean(abs(RH_B-LH_B))./( mean(abs(RH_B).^2)+mean(abs(LH_B).^2) ).^.5;

return