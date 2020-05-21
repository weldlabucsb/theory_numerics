function [xtob,q,Ex] = posToBand1D(V,x)
% POSTOBAND1D Computes projection from position basis to energy bands
%   [xtob,q,Ex] = posToBand1D(V,x) diagonalizes the finite difference
%   representation of the cosine lattice Hamiltonian with lattice depth V
%   over the position vector x. xtob is the unitary transformation into the
%   diagonalized basis. The subspace is size N/2 if N is the length of x. q
%   indexes the states and Ex are the energies.

N = length(x);
dx = x(2)-x(1); D = max(x)-min(x);
Tmat = -gallery('tridiag',N,1,-2,1)/dx^2;
Vmat = V/2*sparse(1:N,1:N,cos(2*x));
[btox,E] = eigs(Tmat+Vmat,N/2,'smallestreal');
E = diag(E);
xtob = btox';
Ex = E(1:N/2);
k = 2*pi/D*(-N/4:N/4-1)';
q = sort(k,'comparisonmethod','abs');

end