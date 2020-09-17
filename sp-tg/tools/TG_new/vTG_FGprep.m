function FG=vTG_FGprep(z,u,v,n2,Av,Ah,Kv,Kh)
% Stability analysis for viscous, diffusive, stratified, parallel shear flow using Fourier-Galerkin method
% Calculate the items needed by "vTG_FG" previously
% Inputs:
% z:    vertical coordinate vector
% u:	velocity component in x direction
% v:	velocity component in y direction
% n2:   stratification
% Av:   vertical eddy viscosity
% Ah:   horizontal eddy viscosity
% Kv:   vertical eddy diffusivity
% Kh:   horizontal eddy diffusivity
% Outputs:
% FG:   structure containing various Fourier integral matrices for use in vTG_FG.
%-------------------------------------------------------------------------%
% Qiang Lian, Xiamen University, China
% Bill Smyth, Oregon State University, USA
% Zhiyu Liu, Xiamen University, Chian
%-------------------------------------------------------------------------%

z=z-z(1);
H=z(end)-z(1);
nz=length(z);
N=floor(nz/2);
D = [1:N]'*pi/H;

F = sqrt(2/H)*sin(D*z')';
G = sqrt(2/H)*cos(D*z')';

for n = 1:N
    FuF(:,n) = trapz(z,bsxfun(@times,bsxfun(@times,F,u),F(:,n)))';
    FvF(:,n) = trapz(z,bsxfun(@times,bsxfun(@times,F,v),F(:,n)))';
	FuzzF(:,n) = ((2*D(n)*D).*trapz(z,bsxfun(@times,bsxfun(@times,G,u),G(:,n)))'-(D(n)^2+D.^2).*trapz(z,bsxfun(@times,bsxfun(@times,F,u),F(:,n)))');
    FvzzF(:,n) = ((2*D(n)*D).*trapz(z,bsxfun(@times,bsxfun(@times,G,v),G(:,n)))'-(D(n)^2+D.^2).*trapz(z,bsxfun(@times,bsxfun(@times,F,v),F(:,n)))');
	FAvzzF(:,n) = ((2*D(n)*D).*trapz(z,bsxfun(@times,bsxfun(@times,G,Av),G(:,n)))'-(D(n)^2+D.^2).*trapz(z,bsxfun(@times,bsxfun(@times,F,Av),F(:,n)))');
	FAvzG(:,n) = (D(n)*trapz(z,bsxfun(@times,bsxfun(@times,F,Av),F(:,n)))'-D.*trapz(z,bsxfun(@times,bsxfun(@times,G,Av),G(:,n)))');
	FAhzAvzG(:,n) = (D(n)*trapz(z,bsxfun(@times,bsxfun(@times,F,Ah+Av),F(:,n)))'-D.*trapz(z,bsxfun(@times,bsxfun(@times,G,Ah+Av),G(:,n)))');
	FAvF(:,n) = trapz(z,bsxfun(@times,bsxfun(@times,F,Av),F(:,n)))';
	FAhAvF(:,n) = trapz(z,bsxfun(@times,bsxfun(@times,F,Ah+Av),F(:,n)))';
	FAhF(:,n) = trapz(z,bsxfun(@times,bsxfun(@times,F,Ah),F(:,n)))';
	FKvzG(:,n) = (D(n)*trapz(z,bsxfun(@times,bsxfun(@times,F,Kv),F(:,n)))'-D.*trapz(z,bsxfun(@times,bsxfun(@times,G,Kv),G(:,n)))');
	FKvF(:,n) = trapz(z,bsxfun(@times,bsxfun(@times,F,Kv),F(:,n)))';
	FKhF(:,n) = trapz(z,bsxfun(@times,bsxfun(@times,F,Kh),F(:,n)))';
	FBzF(:,n) = trapz(z,bsxfun(@times,bsxfun(@times,F,n2),F(:,n)))';
end
% outArgs = {D,F,G,FuF,FvF,FuzzF,FvzzF,FAvzzF,FAvzG,FAhzAvzG,FAvF,FAhAvF,FAhF,FKvzG,FKvF,FKhF,FBzF};
% [varargout{1:nargout}] = outArgs{1:nargout};

FG.D=D;FG.F=F;FG.G=G;FG.FuF=FuF;FG.FvF=FvF;
FG.FuzzF=FuzzF;FG.FvzzF=FvzzF;FG.FAvzzF=FAvzzF;FG.FAvzG=FAvzG;FG.FAhzAvzG=FAhzAvzG;
FG.FAvF=FAvF;FG.FAhAvF=FAhAvF;FG.FAhF=FAhF;FG.FKvzG=FKvzG;FG.FKvF=FKvF;FG.FKhF=FKhF;FG.FBzF=FBzF;


end