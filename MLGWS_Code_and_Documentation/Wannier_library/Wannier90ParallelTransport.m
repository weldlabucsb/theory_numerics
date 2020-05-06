
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
% 1. N. Marzari and D. Vanderbilt. Maximally localized generalized Wannier functions for
% composite energy bands. Phys. Rev. B 56(20), 12847-12865 Nov 1997.
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch. 
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

function [U, omega, omegaI, omegaD, omegaOD, centre, Mmn] = Wannier90ParallelTransport(Mmn, neighbours, mesh)

numBands = size(Mmn, 1); numQpts = size(neighbours.Nearest, 1); numBvecs = length(neighbours.Weight);
dim = size(neighbours.B, 1);
U = eye(numBands); U = U(:, :, ones(1, numQpts)); U1 = U;

% Minimise omegaOD using the parallel transport method, Sec IV.C.1 of Ref. 1
for q = 1 : numQpts % Calculate the transformation matrices
    [v, sigma, w] = svd(Mmn(:,:,q,1));
    U1(:, :, neighbours.Nearest(q, 1)) = w * v';
    U = multiprod(U, U1);
    Mmn = multiprod(conj(permute(U1(:,:,:,ones(1, numBvecs)), [2 1 3 4])), Mmn);
    Mmn = multiprod(Mmn, reshape(U1(:,:,neighbours.Nearest), numBands, numBands, numQpts, numBvecs));
    U1 = eye(numBands); U1 = U1(:,:,ones(1, numQpts));
end
[V2, sigma2] = eig(U(:,:,1));
U1 = V2(:,:,ones(1, numQpts));
U = multiprod(U, U1);
Mmn = multiprod(conj(permute(U1(:,:,:,ones(1, numBvecs)), [2 1 3 4])), Mmn);
Mmn = multiprod(Mmn, reshape(U1(:,:,neighbours.Nearest), numBands, numBands, numQpts, numBvecs));

% Adjust the phases to set omegaD = 0
Mnn = zeros(numBands, numBands, numQpts, numBvecs); % Overlaps between nn. mesh points (within bands)
for n = 1 : numBands
    Mnn(1,n,:,:) = Mmn(n,n,:,:); end
Mnn = Mnn(ones(1, numBands),:,:,:); % n,n,q,b
phase = imag(log(permute(Mnn, [3 4 2 1]))); % Phases between nn. mesh points (within bands), q,b,n,n
phase = permute(phase(:,:,:,1), [3 2 1]); % n,b,q
U2 = ones(numBands, numQpts);
for n = 1 : numBands
    for q = 1 : numQpts - 1 % Calculate the transformation matrices
        U2(n, q+1) = U2(n,q) * exp(-1i * phase(n,1,q)) * exp(1i * sum(phase(n,1,:)) / numQpts); end
end
for q = 1 : numQpts
    U1(:,:,q) = diag(U2(:,q)); end
% Update the transformation matrices and overlap matrices
U = multiprod(U, U1);
Mmn = multiprod(conj(permute(U1(:,:,:,ones(1, numBvecs)), [2 1 3 4])), Mmn);
Mmn = multiprod(Mmn, reshape(U1(:,:,neighbours.Nearest), numBands, numBands, numQpts, numBvecs));

% Calculate the diagonal part of the total spread for the composite group, omegaD
Mnn = zeros(numBands, numBands, numQpts, numBvecs); % Overlaps between nn. mesh points (within bands)
for n = 1 : numBands
    Mnn(1,n,:,:) = Mmn(n,n,:,:); end
Mnn = Mnn(ones(1, numBands),:,:,:); % n,n,q,b
phase = imag(log(permute(Mnn, [3 4 2 1]))); % Phases between nn. mesh points (within bands), q,b,n,n
phase = permute(phase(:,:,:,1), [3 2 1]); % n,b,q
wbphase = neighbours.Weight(ones(1, numBands),:) .* sum(phase, 3); % Weight the phases in each direction, n,b
% Calculate the position of the Wannier centre, Eq. 31
centre = permute(wbphase(:,:,ones(1, dim)), [3,2,1]) .* neighbours.B(:,:,ones(1, numBands)); % dim,b,n
centre = -sum(permute(centre, [1 3 2]), 3) / numQpts; % dim,n
Brdash = transpose(centre) * neighbours.B; % n,b
qn = phase + Brdash(:,:,ones(1, numQpts)); % Eq. 47, n,b,q
qn = permute(qn, [1 3 2]); % n,q,b
% Calculate omegaD, Eq. 36
omegaD = neighbours.Weight * sum(sum(permute(qn .^ 2, [3 1 2]), 3), 2) / numQpts;

% Calculate the off-diagonal part of the total spread for the composite group, omegaOD
Mmn2 = abs(permute(Mmn, [4 1 2 3])) .^ 2; % Absolute value squared of the overlaps
Mnn2 = abs(permute(Mnn(1,:,:,:), [4 2 3 1])) .^ 2; % Absolute value squared of the overlaps (within bands)
omegaOD = neighbours.Weight * (sum(sum(sum(Mmn2, 4), 3), 2) - sum(sum(Mnn2, 3), 2)) / numQpts; % Eq. 35

% Calculate the invariant part of the total spread for the composite group
Mmn = permute(Mmn, [4,1,2,3]);
omegaI = neighbours.Weight * (numBands - sum(sum(sum(conj(Mmn) .* Mmn, 4), 3), 2) / numQpts); % Eq. 34
Mmn = permute(Mmn, [2,3,4,1]);

% Calculate the final total spread, Wannier centres, and transformation matrices
omega = omegaOD + omegaD + omegaI; % Eq. 13
% Move the Wannier centre to the home cell
Ur = exp(1i * transpose(mesh) * round(2*centre)/2); % q,n
Ur = permute(Ur(:,:,ones(1, numBands)), [3,2,1]);
U = U .* Ur; centre = centre - round(2*centre)/2;
