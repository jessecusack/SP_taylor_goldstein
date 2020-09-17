function [DK]=BaryL(z,K,O)
%  Compute the derivative matrices DK by using Barycentric
%  Lagrange interpolation, local interpolated
%  May 01, 2018
%
%  Input:
%  z:        Vertical coordinate vector (can be equal or unequal space)
%  K:        Kth derivative
%  O:        Oth order
%-------------------------------------------------------------------------%
% Qiang Lian, Xiamen University, China
% Bill Smyth, Oregon State University, USA
% Zhiyu Liu, Xiamen University, Chian
%-------------------------------------------------------------------------%
N = length(z);
num = K+O;          % stencil width (odd or even)
if mod(num,2)~=0    % (odd)
    a = fix(num/2);
    j1 = 0;
    for i = 1:a
        j1 = j1+1;
        indx = 1:num;
        zc = z(indx);
        [vdZ,w] = weight(zc);
        D = zeros(1,num); D(j1) = 1;
        for k = 1:K
            D = k*vdZ(j1,:).*((w'/w(j1))*D(j1)-D);
            D(j1) = -sum(D);
            DK(i,indx) = D;
        end 
    end
    
    for i = a+1:N-a
        indx = i-a:i+a;
        zc = z(indx);
        [vdZ,w] = weight(zc);
        D = zeros(1,length(zc)); D(a+1) = 1;
        for k = 1:K
            D = k*vdZ(a+1,:).*((w'/w(a+1))*D(a+1)-D);
            D(a+1) = -sum(D);
            DK(i,indx) = D;
        end 
    end
    
    j3 = num-a;
    for i = N-a+1:N
        j3 = j3+1;
        indx = N-num+1:N;
        zc = z(indx);
        [vdZ,w] = weight(zc);
        D = zeros(1,length(zc)); D(j3) = 1;
        for k = 1:K
            D = k*vdZ(j3,:).*((w'/w(j3))*D(j3)-D);
            D(j3) = -sum(D);
            DK(i,indx) = D;
        end 
    end
else                    %(even)
    a = fix((num-1)/2);
    for i = 1:a
        indx = 1:num;
        zc = z(indx);
        [vdZ,w] = weight(zc);
        D = zeros(1,num); D(i) = 1;
        for k = 1:K
            D = k*vdZ(i,:).*((w'/w(i))*D(i)-D);
            D(i) = -sum(D);
            DK(i,indx) = D;
        end 
    end
    
    for i = a+1:N-a
        nump = num-1;
        indx = i-a:i+a;
        zc = z(indx);
        [vdZ,w] = weight(zc);
        D = zeros(1,nump); D(a+1) = 1;
        for k = 1:K
            D = k*vdZ(a+1,:).*((w'/w(a+1))*D(a+1)-D);
            D(a+1) = -sum(D);
            DK(i,indx) = D;
        end         
    end
    
    j3 = num-a;
    for i = N-a+1:N
        j3 = j3+1;
        indx = N-num+1:N;
        zc = z(indx);
        [vdZ,w] = weight(zc);
        D = zeros(1,num); D(j3) = 1;
        for k = 1:K
            D = k*vdZ(j3,:).*((w'/w(j3))*D(j3)-D);
            D(j3) = -sum(D);
            DK(i,indx) = D;
        end
    end
end

return
end


function [vdZ,w] = weight(z)
    N = length(z);
    I = eye(N);                          % Identity matrix.     
    L = logical(I);                      % Logical identity matrix.

    Z=repmat(z,1,N);
    dZ=Z-Z';
    vdZ = 1./dZ;                         % dZ contains entries 1/(z(i)-z(j))  
    vdZ(L) = zeros(N,1);                 % with zeros on the diagonal.

    % Compute the weights
    w=1./prod(dZ+eye(N),2);


return
end