function [E,ck] = bloch1D(V,q,n)
% BLOCH1D Calculate Bloch bands and Bloch states
%   E = bloch1D(V,q,n) Computes the first n energy bands for cosine lattice
%   of depth V at quasimomentum q. If q is a scalar, E is a vector of
%   length n. If q is a vector, E is an n by length(q) array, where each
%   row represents the energy band over the range q.
%
%   [E,ck] = bloch1D(V,q,n) Also returns the Bloch states in the plane wave
%   basis ck. If q is a scalar, ck is an n by n matrix where the columns
%   correspond to the Bloch states for different band index. If q is a
%   vector, ck is an n by n by length(q) array where the last dimension
%   indexes the quasimomentum.

j = 1-n:2:n-1;
Vmat = V/4*gallery('tridiag',n,1,0,1);
switch isscalar(q)
    case 1
        Tmat = sparse(1:n,1:n,(q+j).^2,n,n);
        switch nargout
            case 1
                E = eig(Vmat+Tmat);
            case 2
                [ck,E] = eig(full(Vmat+Tmat));
                E = diag(E);
        end
    case 0
        E = zeros(n,length(q));
        switch nargout
            case 1
                for ii = 1:length(q)
                    Tmat = sparse(1:n,1:n,(q(ii)+j).^2,n,n);
                    E(:,ii) = eig(Vmat+Tmat);
                end
            case 2
                ck = zeros(n,n,length(q));
                for ii = 1:length(q)
                    Tmat = sparse(1:n,1:n,(q(ii)+j).^2,n,n);
                    [ck(:,:,ii),tempE] = eig(full(Vmat+Tmat));
                    E(:,ii) = diag(tempE);
                end
        end
end
end