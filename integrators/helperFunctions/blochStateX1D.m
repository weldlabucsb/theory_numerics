function unk = blochStateX1D(V,n,q,x)
% BLOCHSTATEX1D Returns the Bloch state u_{n,q}(x)
% DEPENDENCIES: BLOCH1D
%   unk = blochStateX1D(V,n,q,x) computes the Bloch state associated with
%   quasimomentum q and band n for a cosine lattice of depth V. x is an
%   optional position vector; if given, unk is the Bloch wavefunction
%   evaluted over x. If not, unk is an anonymous function handle
%   corresponding to u_{n,q}(x).

nmax = 2*n+21;
[~,v] = bloch1D(V,q,nmax);
k = 1-nmax:2:nmax-1;
vn = v(:,n);
u = @(x) 0;
for ii = 1:nmax
    u = @(x) u(x) + vn(ii)*exp(1i*(k(ii)+q)*x);
end

switch nargin
    case 3
        unk = @(x) u(x);
    case 4
        unk = u(x);
end

end