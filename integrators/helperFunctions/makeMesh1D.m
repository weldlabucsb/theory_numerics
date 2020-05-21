function mesh = makeMesh1D(N,varargin)
% MAKEMESH1D Create 1D space-momentum mesh grid. Must specify the number of
% mesh points N to be an even number, and either the spacing 'delta' or the
% half width 'bound'. To make a grid asymmetric about 0, add a 'shift'.
% Default 'shift' is 0. Returns a struct with the fields:
%   x,k - position and conjugate momentum column vectors.
%   dx,dk - position and momentum spacing
%   Examples of valid calls:
%   mesh = makeMesh1D(N,'delta',dx) - from -dx*N/2 to dx*N/2-1
%   mesh = makeMesh1D(N,'bound',bound,'shift',shift) - from -bound+shift to
%   bound+shift

    p = inputParser;
    isPosScalar = @(x) isnumeric(x) && isscalar(x) && (x>0);
    checkN = @(x) isPosScalar(x) && mod(x,2)==0;
    addRequired(p,'N',checkN);
    addOptional(p,'delta',0,isPosScalar);
    addOptional(p,'bound',0,isPosScalar);
    addOptional(p,'shift',0,@(x) isnumeric(x) && isscalar(x));
    parse(p,N,varargin{:});
    
    if floor(log2(p.Results.N)) ~= ceil(log2(p.Results.N))
        warning('Number of mesh points not a power of 2.')
    end
    
    switch p.Results.delta
        case 0
            if p.Results.bound == 0
                error('Mesh size 0. Not enough input arguments specified.')
            end
            x = linspace(-p.Results.bound,p.Results.bound,N)'+p.Results.shift;
            dx = x(2)-x(1);
            dk = pi/p.Results.bound;
            k = dk*(-p.Results.N/2:p.Results.N/2-1)';
        otherwise
            if p.Results.bound ~= 0
                warning('Mesh overspecified. Creating mesh based on delta.')
            end
            dx = p.Results.delta;
            x = p.Results.shift+dx*(-p.Results.N/2:p.Results.N/2-1)';
            dk = 2*pi/(max(x)-min(x));
            k = dk*(-p.Results.N/2:p.Results.N/2-1)';
    end
    
    mesh.x = x;
    mesh.k = k;
    mesh.dx = dx;
    mesh.dk = dk;
end