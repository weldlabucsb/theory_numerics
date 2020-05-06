
% Calculate maximally localised Wannier states for an isolated band (lattice with inversion symmetry)
% 
% U: transformation matrices acting on the Bloch states to produce maximally localised Wannier
%    functions.
% Mmn: updated overlaps between neighbouring q-mesh points.
% omega: Spread functional of the Wannier function.
% omegaI: invariant part of the spread functional.
% omegaD: diagonal part of the spread functional.
% rdash: mean position of the Wannier functions.
% 
% REFERENCES (equation numbers are referenced in the code):
% 1. N. Marzari and D. Vanderbilt. Maximally localized generalized Wannier functions for
% composite energy bands. Phys. Rev. B 56(20), 12847-12865 Nov 1997.
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

function [U, Mmn, omega, omegaI, omegaD, centre] = Wannier90Inversion(dimension, neighbours, Mmn)

numQpts = size(neighbours.Nearest, 1); numBvecs = length(neighbours.Weight); dim = size(neighbours.B, 1);
U = ones(numQpts, 1);

% Calculate the phases between nearest-neighbour q-points
phase = imag(log(Mmn));

% Sweep through the unit cell, equalising the phases between neighbours
switch dimension
    case '1D'
        for q = 1 : numQpts - 1 % Calculate the transformation matrices
            U(q + 1) = U(q) * exp(-1i * phase(q, 1)) * exp(1i * sum(phase(:, 1)) / numQpts); end
    case '2D'
        N = sqrt(numQpts);
        % Calculate the Berry phase along the direction of the first reciprocal lattice vector
        phasex = reshape(phase(:,1), N, N);
        sumphasex = sum(phasex, 2) / (2 * pi);
        sumphasex = abs(2 * pi * (sumphasex - round(sumphasex)));
        % Calculate the Berry phase along the direction of the second reciprocal lattice vector
        phasey = reshape(phase(:,3), N, N);
        sumphasey = sum(phasey, 1) / (2 * pi);
        sumphasey = abs(2 * pi * (sumphasey - round(sumphasey)));
        q = 1; % Calculate the transformation matrices
        for q1 = 1 : N - 1
            for q2 = 1 : N - 1
                U(q + 1) = U(q) * exp(-1i * phasey(q2, q1)) * exp(1i * sumphasey(q1) / N); q = q + 1; end
            U(q + 1) = U(q + 1 - N) * exp(-1i * phasex(1, q1)) * exp(1i * sumphasex(1) / N); q = q + 1;
        end
        for q2 = 1 : N - 1
            U(q + 1) = U(q) * exp(-1i * phasey(q2, N)) * exp(1i * sumphasey(N) / N); q = q + 1; end
    case '3D'
        N = round(numQpts ^ (1 / 3));
        % Calculate the Berry phase along the direction of the first reciprocal lattice vector
        phasex = reshape(phase(:,1), N, N, N);
        sumphasex = permute(sum(phasex, 2) / (2 * pi), [1 3 2]);
        sumphasex = abs(2 * pi * (sumphasex - round(sumphasex)));
        % Calculate the Berry phase along the direction of the second reciprocal lattice vector
        phasey = reshape(phase(:,3), N, N, N);
        sumphasey = permute(sum(phasey, 1) / (2 * pi), [2 3 1]);
        sumphasey = abs(2 * pi * (sumphasey - round(sumphasey)));
        % Calculate the Berry phase along the direction of the third reciprocal lattice vector
        phasez = reshape(phase(:,5), N, N, N);
        sumphasez = sum(phasez, 3) / (2 * pi);
        sumphasez = abs(2 * pi * (sumphasez - round(sumphasez)));
        q = 1; % Calculate the transformation matrices
        for q3 = 1 : N - 1
            for q1 = 1 : N - 1
                for q2 = 1 : N - 1
                    U(q + 1) = U(q) * exp(-1i * phasey(q2,q1,q3)) * exp(1i * sumphasey(q1,q3) / N);
                    q = q + 1;
                end
                U(q + 1) = U(q + 1 - N) * exp(-1i * phasex(1,q1,q3)) * exp(1i * sumphasex(1,q3) / N);
                q = q + 1;
            end
            for q2 = 1 : N - 1
                U(q + 1) = U(q) * exp(-1i * phasey(q2,N,q3)) * exp(1i * sumphasey(N,q3) / N); q = q + 1; end
            U(q + 1) = U(q + 1 - N ^ 2) * exp(-1i * phasez(1,1,q3)) * exp(1i * sumphasez(1,1) / N); q = q + 1;
        end
        for q1 = 1 : N - 1
            for q2 = 1 : N - 1
                U(q + 1) = U(q) * exp(-1i * phasey(q2,q1,N)) * exp(1i * sumphasey(q1,N) / N); q = q + 1; end
            U(q + 1) = U(q + 1 - N)*exp(-1i * phasex(1,q1,N)) * exp(1i * sumphasex(1,N) / N); q = q + 1;
        end
        for q2 = 1 : N - 1
            U(q + 1) = U(q) * exp(-1i * phasey(q2,N,N)) * exp(1i * sumphasey(N,N) / N); q = q + 1; end
end

% Realculate the overlaps
Mmn = conj(U(:, ones(1, numBvecs))) .* Mmn .* U(neighbours.Nearest); % Eq. 61

% Calculate the diagonal part of the total spread, omegaD
phase = imag(log(Mmn)); % The phases between nearest-neighbour mesh points
% Calculate the position of the Wannier centre, Eq. 31
centre = -sum(repmat(neighbours.Weight .* sum(phase, 1), dim, 1) .* neighbours.B, 2) / numQpts;
Brdash = transpose(centre) * neighbours.B; qn = phase + Brdash(ones(1, numQpts), :); % Eq. 47
omegaD = sum(((-qn) .^ 2) * transpose(neighbours.Weight)) / numQpts; % Calculate omegaD, Eq. 36

% Calculate the invariant part of the total spread
omegaI = neighbours.Weight * transpose(1 - sum(abs(Mmn) .^ 2, 1) / numQpts); % Eq. 34

% Calculate the final total spread, Wannier centres, and transformation matrices
omega = omegaD + omegaI; % Eq. 13

