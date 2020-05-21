function c = strangNI(c,Uk,UV1,UV2)
% STRANGNI Performs Strang propagation over a single time step for
% noninteracting system.
%   c = strangNI(c,Uk,UV1,UV2) returns the evolution of the wavefunction defined
%   by c over a small time step, given the momentum propagator Uk over the
%   full step and position propagators UV1 and UV2 over the half steps.

if ~isequal(size(c),size(Uk),size(UV1),size(UV2))
    error('Inconsistent input dimensions');
end

switch ndims(c)
    case 2
        switch isvector(c)
            case 1
                c = UV1.*c;
                c = ifft(ifftshift(Uk.*fftshift(fft(c))));
                c = UV2.*c;
            case 0
                c = UV1.*c;
                c = ifft2(ifftshift(Uk.*fftshift(fft2(c))));
                c = UV2.*c;
        end
    otherwise
        c = UV1.*c;
        c = ifftn(ifftshift(Uk.*fftshift(fftn(c))));
        c = UV2.*c;
end
end