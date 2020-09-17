function varargout=vTG_FG(z,u,v,Av,Ah,k,l,imode,FG)
% Stability analysis for viscous, diffusive, stratified, parallel shear flow using Fourier-Galerkin method
% Inputs:
% z:    vertical coordinate vector
% u:	velocity component in x direction
% v:	velocity component in y direction
% k:    wavenumber in x direction
% l:    wavenumber in y direction
% imode:mode selection in terms of growth rate
%       0 = all modes
%       1 = fastest-growing mode
% FG:   structure containing Fourier variables computed previously in vTG_FGprep.
%
%-------------------------------------------------------------------------%
% Qiang Lian, Xiamen University, China
% Bill Smyth, Oregon State University, USA
% Zhiyu Liu, Xiamen University, Chian
%-------------------------------------------------------------------------%
ii = complex(0,1);
z0 = z;
z=z-z(1);
nz=length(z);
N=floor(nz/2);

D=FG.D;
F=FG.F;
G=FG.G;
FuF=FG.FuF;
FvF=FG.FvF;
FuzzF=FG.FuzzF;
FvzzF=FG.FvzzF;
FAvzzF=FG.FAvzzF;
FAvzG=FG.FAvzG;
FAhzAvzG=FG.FAhzAvzG;
FAvF=FG.FAvF;
FAhAvF=FG.FAhAvF;
FAhF=FG.FAhF;
FKvzG=FG.FKvzG;
FKvF=FG.FKvF;
FKhF=FG.FKhF;
FBzF=FG.FBzF;

D2 = D.^2; 
D3 = D.^3; 
D4 = D.^4;

k2 = k^2+l^2;
kt=sqrt(k2);

U = u*(k/kt)+v*(l/kt);    % velocity component parallel to the wave vector (k,l)  
FUF = FuF*(k/kt)+FvF*(l/kt);
FUzzF = FuzzF*(k/kt)+FvzzF*(l/kt);

lap=-D2-k2;
B11 = (-ii*kt*bsxfun(@times,FUF,lap')+(ii*kt*FUzzF)-bsxfun(@times,FAvzzF,D2')-...
        2*bsxfun(@times,FAvzG,D3')+bsxfun(@times,FAvF,D4')-k2*bsxfun(@times,FAhzAvzG,D')+...
        k2*bsxfun(@times,FAhAvF,D2')+k2^2*FAhF)./(lap*ones(1,N));
B21 =  FBzF;
B22 = ((-ii*kt*FUF)+bsxfun(@times,FKvzG,D')-bsxfun(@times,FKvF,D2')-(k2*FKhF));
RHS=[B11 -diag(k2./lap);-B21 B22];
        
[X,E]=eig(RHS);
E=diag(E);
%rem = abs(E)>1000; E(rem) = []; X(:,rem) = [];
        
% Sort eigvals and eigvecs by real growth rate
[~,ind]=sort(real(E),1,'descend');
X=X(:,ind);
E=E(ind);

% Save the mode selected via imode
if imode==0
	sig=E;
	w=F*X(1:N,:);
	b=F*X(N+1:end,:);
elseif imode>0
	sig=E(imode);
	wf=X(1:N,imode); bf=X(N+1:end,imode);
	w=F*wf; b=F*bf;
    
    % normalize
 	cnorm=w(find(abs(w)==max(abs(w)),1));
 	w=w/cnorm;wf=wf/cnorm; 
 	b=b/cnorm;
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	Kinetic energy budget

    wz=(diag(D)*G')'*wf; 
    wzz = -(diag(D2)*F')'*wf; 
    wzzz = -(diag(D3)*G')'*wf;
	u=(ii/kt)*wz;
	uz=(ii/kt)*wzz;
	uzz=(ii/kt)*wzzz;
	K=(abs(u).^2+abs(w).^2)/4;
    
    Bary=BaryL(z0,1,6);
    Uz=Bary*U;
    Avz=Bary*Av;
        
% 	Compute shear production terms
	uw=real(u.*conj(w))/2;
	SP=-uw.*Uz;

% 	Compute pressure eigfn and energy flux
	p=-(sig+ii*kt*U+Ah*k2).*wz/k2 +(ii/kt)*Uz.*w+(Av/k2).*(wzzz)+(Avz/k2).*(wzz);
	pw=real(p.*conj(w))/2;
	EFz=Bary*pw;

%   Compute buoyancy flux
	BF=real(b.*conj(w))/2;
    
% 	Compute viscous effects
    AvKz_z = .5*Avz.*real(conj(u).*uz+conj(w).*wz) ...
            +.5*Av.*(abs(uz).^2+abs(wz).^2+real(conj(u).*uzz+conj(w).*wzz));
	Eps=0.5*Av.*(abs(uz).^2+abs(wz).^2)+2*Ah.*k2.*K; 
end

% The output:
if ( nargout <= 1 && imode==0)
	varargout{1} = sig;
elseif ( nargout > 1 && imode==0)
    outArgs = {sig,w,b};
	[varargout{1:nargout}] = outArgs{1:nargout};
elseif ( nargout <= 13 && imode>0 )
	outArgs = {sig,w,b,uw,SP,pw,EFz,BF,AvKz_z,K,p,Eps};
	[varargout{1:nargout}] = outArgs{1:nargout};
else
	error('vTG_FG:nargout', ...
        'Incorrect number of output arguments.'); 
end 

end