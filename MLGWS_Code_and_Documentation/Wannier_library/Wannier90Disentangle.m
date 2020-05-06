
% Disentangle a group of bands sequentially by minimising omegaI for each band in turn
% 
% U: transformation matrices acting on the Bloch states to disentangle the bands.
% Mmn: updated overlaps between neighbouring q-mesh points.
% 
% REFERENCES (equation numbers are referenced in the code):
% 1. I. Souza, N. Marzari and D. Vanderbilt. Maximally localized Wannier functions for
% entangled energy bands. Phys. Rev. B 65(3), 035109 Dec 2001.
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch. 
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

function [U, omega, omegaI, omegaD, omegaOD, centre, Mmn] = Wannier90Disentangle(Mmn, neighbours, ...
    alpha, maxIter, randomise)

numQpts = size(neighbours.Nearest, 1); numBvecs = length(neighbours.Weight); numBands = size(Mmn, 1);
dim = size(neighbours.B, 1);
U = eye(numBands); U = U(:,:,ones(1, numQpts));

% Randomise the band labels at each mesh point
if strcmp(randomise, 'true')
    for q = 1 : numQpts
        U1 = eye(numBands); U1 = U1(randperm(numBands), :); U(:,:,q) = U(:,:,q) * U1; end
    Mmn = multiprod(conj(permute(U(:,:,:,ones(1, numBvecs)), [2 1 3 4])), Mmn);
    Mmn = multiprod(Mmn, reshape(U(:,:,neighbours.Nearest), numBands, numBands, numQpts, numBvecs));
end

Mmni1 = Mmn; % For the first iteration, set the input overlaps for the previous loop equal to the
% input overaps for the current loop 
for n = 1 : numBands - 1
    UTemp = zeros(numBands - n + 1, numBands - n + 1, numQpts); eigval = zeros(numBands - n + 1, numQpts);
    Mmn1 = permute(Mmn, [3,4,1,2]); loop = 1; omegaIn = zeros(1, maxIter + 1);
    omegaIn(loop) = neighbours.Weight * transpose(1 - sum(abs(Mmn1(:,:,n,n)) .^ 2, 1) / numQpts);
    disp(['n = ' num2str(n) '/' num2str(numBands - 1) ': iteration = 0, omegaI = ' num2str(omegaIn(loop), 6)]);
    for loop = 1 : maxIter
        U1 = eye(numBands); U1 = U1(:,:,ones(1, numQpts));
        % Calculate the Z matrix, Eq. 21
        Mmn2 = permute(Mmn(n:numBands,n,:,:), [1 4 3 2]) .* neighbours.Weight(ones(1, numBands - n + 1), ...
            :, ones(1, numQpts));
        Mmn2 = multiprod(Mmn2, conj(permute(Mmn(n:numBands,n,:,:), [4 1 3 2])));
        Mmn3 = permute(Mmni1(n:numBands,n,:,:), [1 4 3 2]) .* neighbours.Weight(ones(1, numBands - n + 1), ...
            :, ones(1, numQpts));
        Mmn3 = multiprod(Mmn3, conj(permute(Mmni1(n:numBands,n,:,:), [4 1 3 2])));
        Z = alpha * Mmn2 + (1 - alpha) * Mmn3;
        for q = 1 : numQpts % Diagonalise the Z matrix
            [UTemp(:,:,q), eigvalTemp] = eig(Z(:,:,q));
            eigval(:,q) = real(diag(eigvalTemp));
        end
        % Sort the eigenvalues of Z into descending order, Eq. 19
        [eigvalTemp2, eigInd] = sort(eigval,'descend'); 
        eigInd2 = 0 : numBands - n + 1 : (numQpts - 1) * (numBands - n + 1);
        eigInd2 = eigInd2(ones(1, numBands - n + 1 ), :);
        eigInd = eigInd + eigInd2;
        % Sort the corresponding transformation matrices according to the order of the eigenvalues
        UTemp = reshape(UTemp(:, eigInd), numBands - n + 1, numBands - n + 1, numQpts);
        % Update the full set of transformation matrices
        if n == 1
            U1 = multiprod(U1, UTemp);
        else
            U2 = eye(numBands); U2 = U2(:,:,ones(1, numQpts));
            U2(n:numBands, n:numBands, :) = UTemp;
            U1 = multiprod(U1, U2);
        end
        % Update the overlaps between nn. mesh points, Eq. 61
        Mmn = multiprod(conj(permute(U1(:,:,:,ones(1, numBvecs)), [2 1 3 4])), Mmn);
        Mmni1 = Mmn; % Input overlaps for the previous loop, Eq. 20
        Mmn = multiprod(Mmn, reshape(U1(:,:,neighbours.Nearest), numBands, numBands, numQpts, numBvecs));
        U = multiprod(U, U1); % Update the full set of transformation matrices
        Mmn1 = permute(Mmn, [3,4,1,2]);
        % Recalculate the invariant part of the spread, omegaI for the current band, Eq. 7
        omegaIn(loop + 1) = neighbours.Weight * transpose(1 - sum(abs(Mmn1(:,:,n,n)) .^ 2, 1) / numQpts);
        if (loop) / 10 == round((loop) / 10)
            disp(['n = ' num2str(n) '/' num2str(numBands - 1) ': iteration ' num2str(loop) ', omegaI = ' ...
                num2str(omegaIn(loop + 1), 6)]);
        end
    end
end

% -------------------------------------------------------------------------------
% Calculate the diagonal part of the total spread for the composite group, omegaD
% -------------------------------------------------------------------------------
Mnn = zeros(numBands, numBands, numQpts, numBvecs); % Overlaps between nn. mesh points (within bands)
for n = 1 : numBands
    Mnn(1,n,:,:) = Mmn(n,n,:,:); end
Mnn = Mnn(ones(1, numBands),:,:,:); % n,n,q,b
phase = imag(log(permute(Mnn, [3 4 2 1]))); % Phases between nn. mesh points (within bands), q,b,n,n
phase = permute(phase(:,:,:,1), [3 2 1]); % n,b,q
wbphase = neighbours.Weight(ones(1, numBands), :) .* sum(phase, 3); % Weight the phases in each direction, n,b
% Calculate the position of the Wannier centre, Eq. 31
centre = permute(wbphase(:, :, ones(1, dim)), [3,2,1]) .* neighbours.B(:, :, ones(1, numBands)); % dim,b,n
centre = -sum(permute(centre, [1 3 2]), 3) / numQpts; % dim,n
Brdash = transpose(centre) * neighbours.B; % n,b
qn = phase + Brdash(:, :, ones(1, numQpts)); % Eq. 47, n,b,q
qn = permute(qn, [1 3 2]); % n,q,b
% Calculate omegaD, Eq. 36
omegaD = neighbours.Weight * sum(sum(permute(qn .^ 2, [3 1 2]), 3), 2) / numQpts;

% ------------------------------------------------------------------------------------
% Calculate the off-diagonal part of the total spread for the composite group, omegaOD
% ------------------------------------------------------------------------------------
Mmn2 = abs(permute(Mmn, [4 1 2 3])) .^ 2; % Absolute value squared of the overlaps
Mnn2 = abs(permute(Mnn(1,:,:,:), [4 2 3 1])) .^ 2; % Absolute value squared of the overlaps (within bands)
omegaOD = neighbours.Weight * (sum(sum(sum(Mmn2, 4), 3), 2) - sum(sum(Mnn2, 3), 2)) / numQpts; % Eq. 35

% ------------------------------------------------------------------------
% Calculate the invariant part of the total spread for the composite group
% ------------------------------------------------------------------------
Mmn = permute(Mmn, [4,1,2,3]);
omegaI = neighbours.Weight * (numBands - sum(sum(sum(conj(Mmn) .* Mmn, 4), 3), 2) / numQpts); % Eq. 34
Mmn = permute(Mmn, [2,3,4,1]);

omega = omegaI + omegaD + omegaOD;
