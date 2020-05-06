
% Wannier90 steepest-descent algorithm to produce maximally localised Wannier states (composite group)
% 
% U: transformation matrices acting on the Bloch states to produce maximally localised Wannier functions.
% omega: Spread functional of the Wannier function.
% omegaI: invariant part of the spread functional.
% omegaD: diagonal part of the spread functional.
% omegaOD: off-diagonal part of the spread functional.
% rdash: mean position of the Wannier functions.
% Mmn: updated overlaps between neighbouring q-mesh points.
% 
% REFERENCES (equation numbers are referenced in the code):
% 1. N. Marzari and D. Vanderbilt. Maximally localized generalized Wannier functions for composite energy
% bands. Phys. Rev. B 56(20), 12847-12865 Nov 1997.
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch. 
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

function [U, omega, omegaI, omegaD, omegaOD, centre, Mmn] = Wannier90Composite(Mmn, neighbours, ...
    epsilon, maxIter, mesh)

numBands = size(Mmn, 1); numQpts = size(neighbours.Nearest, 1); numBvecs = length(neighbours.Weight); dim = size(neighbours.B, 1);
epsilon = epsilon / (4 * sum(neighbours.Weight)); omegaD = zeros(1, maxIter + 1); omegaOD = zeros(1, maxIter + 1);
U = eye(numBands); U = U(:, :, ones(1, numQpts)); U1 = U;

% -------------------------------------------------------------------------------
% Calculate the diagonal part of the total spread for the composite group, omegaD
% -------------------------------------------------------------------------------
Mnn = zeros(numBands, numBands, numQpts, numBvecs); % Overlaps between nn. mesh points (within bands)
for n = 1 : numBands
    Mnn(1,n,:,:) = Mmn(n,n,:,:); end
Mnn = Mnn(ones(1, numBands),:,:,:); % n,n,q,b
% Phases between nn. mesh points (within bands), q,b,n,n
phase = imag(log(permute(Mnn, [3 4 2 1])));
phase = permute(phase(:,:,:,1), [3 2 1]); % n,b,q
% Weight the phases in each direction, n,b
wbphase = neighbours.Weight(ones(1, numBands), :) .* sum(phase, 3);
% Calculate the position of the Wannier centre, Eq. 31
centre = permute(wbphase(:, :, ones(1, dim)), [3,2,1]) .* neighbours.B(:, :, ones(1, numBands)); % dim,b,n
centre = -sum(permute(centre, [1 3 2]), 3) / numQpts; % dim,n
Brdash = transpose(centre) * neighbours.B; % n,b
qn = phase + Brdash(:, :, ones(1, numQpts)); % Eq. 47, n,b,q
qn = permute(qn, [1 3 2]); % n,q,b
% Calculate omegaD, Eq. 36
loop = 1; omegaD(loop) = neighbours.Weight * sum(sum(permute(qn .^ 2, [3 1 2]), 3), 2) / numQpts;

% ------------------------------------------------------------------------------------
% Calculate the off-diagonal part of the total spread for the composite group, omegaOD
% ------------------------------------------------------------------------------------
Mmn2 = abs(permute(Mmn, [4 1 2 3])) .^ 2; % Absolute value squared of the overlaps
Mnn2 = abs(permute(Mnn(1,:,:,:), [4 2 3 1])) .^ 2; % Absolute value squared of the overlaps (within bands)
omegaOD(loop) = neighbours.Weight * (sum(sum(sum(Mmn2, 4), 3), 2) - sum(sum(Mnn2, 3), 2)) / numQpts; % Eq. 35
disp(['Iteration = 0, omegaOD = ' num2str(omegaOD(loop), 4) ', omegaD = ' num2str(omegaD(loop), 4)]);

% -------------------------------------------------------------------------
% Use the Wannier90 steepest-descent algorithm to reduce omegaOD and omegaD
% -------------------------------------------------------------------------
while loop <= maxIter
    Rmn = Mmn .* conj(Mnn); % Eq. 45
    T = permute(qn(:,:,:,ones(1, numBands)), [4 1 2 3]) .* Mmn ./ Mnn; % Eq. 51
    % Calculate the gradient of the spread functional, Eq. 52
    RTmn = Rmn - conj(permute(Rmn, [2 1 3 4])) + 1i * (T + conj(permute(T, [2 1 3 4])));
    Gr = 2 * sum(permute(neighbours.Weight(ones(1, numQpts), :, ones(1, numBands), ones(1, numBands)), ...
        [3 4 1 2]) .* RTmn, 4);
    % Steepest-descent transformation matrices, Eq. 60
    for q = 1 : numQpts
        U1(:,:,q) = expm(epsilon * Gr(:,:,q)); end
    U = multiprod(U, U1);
    % Update the overlaps between nn. mesh points, Eq. 61
    Mmn = multiprod(conj(permute(U1(:, :, :, ones(1, numBvecs)), [2 1 3 4])), Mmn);
    Mmn = multiprod(Mmn, reshape(U1(:,:,neighbours.Nearest), numBands, numBands, numQpts, numBvecs));
    % Update the overlaps between nn. mesh points (within bands)
    Mnn = zeros(numBands, numBands, numQpts, numBvecs);
    for n = 1 : numBands
        Mnn(1,n,:,:) = Mmn(n,n,:,:); end
    Mnn = Mnn(ones(1, numBands),:,:,:); % n,n,q,b
    phase = imag(log(permute(Mnn, [3 4 2 1]))); % Update the phases between nn. mesh points, q,b,n,n
    phase = permute(phase(:,:,:,1), [3 2 1]); % n,b,q
    wbphase = neighbours.Weight(ones(1, numBands), :) .* sum(phase, 3); % Weight the phases in each direction, n,b
    % Recalculate the Wannier centre, Eq. 31
    centre = permute(wbphase(:, :, ones(1, dim)), [3,2,1]) .* neighbours.B(:, :, ones(1, numBands)); % dim,b,n
    centre = -sum(permute(centre, [1 3 2]), 3) / numQpts; % dim,n
    Brdash = transpose(centre) * neighbours.B; % n,b
    qn = phase + Brdash(:, :, ones(1, numQpts)); % Recalculate qn, Eq. 47, n,b,q
    qn = permute(qn, [1 3 2]); % n,q,b
    % Recalculate omegaD, Eq. 36
    loop = loop + 1;
    omegaD(loop) = neighbours.Weight * sum(sum(permute(qn .^ 2, [3 1 2]), 3), 2) / numQpts; 
    Mmn2 = abs(permute(Mmn, [4 1 2 3])) .^ 2; % Absolute value squared of the overlaps
    Mnn2 = abs(permute(Mnn(1,:,:,:),[4 2 3 1])) .^ 2; % Absolute value squared of the overlaps (within bands)
    % Recalculate omegaOD, Eq. 35
    omegaOD(loop) = neighbours.Weight * (sum(sum(sum(Mmn2, 4), 3), 2) - sum(sum(Mnn2, 3), 2)) / numQpts;
    if (loop - 1) / 10 == round((loop - 1) / 10)
        disp(['Iteration = ' num2str(loop-  1) ', omegaOD = ' num2str(omegaOD(loop), 4) ...
             ', omegaD = ' num2str(omegaD(loop), 4)]);
    end
end
omegaD(omegaD == 0) = []; omegaOD(omegaOD == 0) = []; % Remove cells if loop<iter
disp(['The total no. of wannier90 iterations is ', num2str(loop - 1)]);

% ------------------------------------------------------------------------
% Calculate the invariant part of the total spread for the composite group
% ------------------------------------------------------------------------
Mmn = permute(Mmn, [4,1,2,3]);
omegaI = neighbours.Weight * (numBands - sum(sum(sum(conj(Mmn) .* Mmn, 4), 3), 2) / numQpts); % Eq. 34
Mmn = permute(Mmn, [2,3,4,1]);

% ------------------------------------------------------------------------------
% Calculate the final total spread, Wannier centres, and transformation matrices
% ------------------------------------------------------------------------------
omega = omegaOD + omegaD + omegaI; % Eq. 13

% Move the Wannier centre to the home cell
% Ur = exp(1i * transpose(mesh) * round(2*centre)/2); % q,n
% Ur = permute(Ur(:,:,ones(1, numBands)), [3,2,1]);
% U = U .* Ur; centre = centre - round(2*centre)/2;


% % Move the Wannier centre to the home cell
% Ur = exp(1i * transpose(mesh) * round(centre)); % q,n
% Ur = permute(Ur(:, :, ones(1, numBands)), [3,2,1]);
% U = U .* Ur;
% U1 = eye(numBands); U1 = U1(:, :, ones(1, numQpts));
% U1 = U1 .* Ur;
% % U = multiprod(U, U1);
% centre = centre - round(centre);
% % Update the overlaps between nn. mesh points, Eq. 61
% Mmn = multiprod(conj(permute(U1(:, :, :, ones(1, numBvecs)), [2 1 3 4])), Mmn);
% Mmn = multiprod(Mmn, reshape(U1(:,:,neighbours.Nearest), numBands, numBands, numQpts, numBvecs));

% % Fix for 1D super-lattice
% centreTemp = centre > 0.4;
% Ur2 = exp(1i * transpose(mesh) * centreTemp); % q,n
% Ur2 = permute(Ur2(:, :, ones(1, numBands)), [3,2,1]);
% U = U .* Ur2;
% U1 = eye(numBands); U1 = U1(:, :, ones(1, numQpts));
% U1 = U1 .* Ur2;
% centre = centre - centreTemp;
% % Update the overlaps between nn. mesh points, Eq. 61
% Mmn = multiprod(conj(permute(U1(:, :, :, ones(1, numBvecs)), [2 1 3 4])), Mmn);
% Mmn = multiprod(Mmn, reshape(U1(:,:,neighbours.Nearest), numBands, numBands, numQpts, numBvecs));
