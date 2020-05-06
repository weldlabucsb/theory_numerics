
% Wannier90 steepest-descent algorithm to produce maximally localised Wannier states (isolated band)
% 
% U: transformation matrices acting on the Bloch states to produce maximally localised Wannier
%    functions.
% omega: Spread functional of the Wannier function.
% omegaI: invariant part of the spread functional.
% omegaD: diagonal part of the spread functional.
% rdash: mean position of the Wannier functions.
% Mmn: updated overlaps between neighbouring q-mesh points.
% 
% REFERENCES (equation numbers are referenced in the code):
% 1. N. Marzari and D. Vanderbilt. Maximally localized generalized Wannier functions for
% composite energy bands. Phys. Rev. B 56(20), 12847-12865 Nov 1997.
% 
% --------------------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch. 
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% --------------------------------------------------------------------------------------------------

function [U, omega, omegaI, omegaD, centre, Mmn] = Wannier90Isolated(Mmn, neighbours, epsilon, maxIter, mesh)

numQpts = size(neighbours.Nearest, 1); numBvecs = length(neighbours.Weight); dim = size(neighbours.B, 1);
epsilon = epsilon / (4 * sum(neighbours.Weight));
U = ones(numQpts, 1); omegaD = zeros(1, maxIter + 1);

% Calculate the initial diagonal part of the total spread, omegaD
phase = imag(log(Mmn)); % Phases between nn. mesh points
centre = -sum(repmat(neighbours.Weight .* sum(phase, 1), dim, 1) .* neighbours.B, 2) / numQpts; % Position of the Wannier centre, Eq. 31
Brdash = transpose(centre) * neighbours.B; qn = phase + Brdash(ones(1, numQpts), :);
loop = 1; omegaD(loop) = sum(((-qn) .^ 2) * transpose(neighbours.Weight)) / numQpts; % Eq. 47
disp(['Iteration = 0, omegaD = ' num2str(omegaD(loop), 4)]); % Eq. 36

% Use the Wannier90 steepest-descent algorithm to reduce omegaD
while loop <= maxIter
    Gr = 4i * phase * transpose(neighbours.Weight); % Gradient of the spread functional, Eq. 52 and Eq. 53
    U1 = exp(epsilon * Gr); U = U .* U1; % Steepest-descent transformation matrices, Eq. 60
    Mmn = conj(U1(:, ones(1, numBvecs))) .* Mmn .* U1(neighbours.Nearest); % Update the overlaps between nn. mesh points, Eq. 61
    phase = imag(log(Mmn)); % Update the phases between nn. mesh points
    centre = -sum(repmat(neighbours.Weight .* sum(phase, 1), dim, 1) .* neighbours.B, 2) / numQpts; % Recalculate the Wannier centre, Eq. 31
    Brdash = transpose(centre) * neighbours.B; qn = phase + Brdash(ones(1, numQpts), :); % Recalculate qn, Eq. 47
    loop = loop + 1; omegaD(loop) = sum(((-qn) .^ 2) * transpose(neighbours.Weight)) / numQpts; % Recalculate omegaD, Eq. 36
    if (loop - 1) / 100 == round((loop - 1) / 100)
        disp(['Iteration = ' num2str(loop - 1) ', omegaD = ' num2str(omegaD(loop), 4)]); end
end
omegaD(omegaD == 0) = []; % Remove cells if loop<iter
if isempty(omegaD)
    omegaD = 0; end
disp(['The total no. of wannier90 iterations is ', num2str(loop - 1)]);

% Calculate the invariant part of the total spread, Eq. 34
omegaI = neighbours.Weight * transpose(1 - sum(abs(Mmn) .^ 2, 1) / numQpts);

% Calculate the final total spread, Wannier centres, and transformation matrices
omega = omegaD + omegaI; % Eq. 13
% Move the Wannier centre to the home cell
U = U .* exp(1i * transpose(mesh) * round(2*centre)/2); centre = centre - round(2*centre)/2;

