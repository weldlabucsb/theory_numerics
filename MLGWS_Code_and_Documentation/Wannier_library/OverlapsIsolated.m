
% Calculate the overlap matrices Mmn for an isolated band
% 
% Mmn: overlaps between neighbouring q-mesh points, labelled by q and b.
% 
% REFERENCES:
% 1. N. Marzari and D. Vanderbilt. Maximally localized generalized Wannier functions for
% composite energy bands. Phys. Rev. B 56(20), 12847-12865 Nov 1997.
%
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

function [Mmn] = OverlapsIsolated(state, meshInd, recip, neighbours)

% ------------------------------------------
% Calculate the overlap matrices Mmn, Eq. 25
% ------------------------------------------
numKpts = size(state, 1); numQpts = size(state, 2); numBvecs = size(neighbours.Nearest, 2);
Mmn = conj(state(:, :, ones(1, numBvecs))) .* reshape(state(:, neighbours.Nearest), numKpts, numQpts, numBvecs);
Mmn = sum(permute(Mmn, [2 3 1]), 3);

% -------------------------------------------------------------------------------------------------------
% Corrections for wrapping around the cell (need to translate sides/edges/corners by a reciprocal lattice
% vector)
% -------------------------------------------------------------------------------------------------------
switch recip.Dimension
    case '1D'
        b = 1; % 1st nearest neighbour vector, G1
        q = numQpts; % Correction for right-hand side of cell
        Mmn(q,b) = state(1 : numKpts - 1, q)' * state(2 : numKpts, neighbours.Nearest(q,b));
        b = 2; % 2nd nearest neighbour vector, -G1
        q = 1; % Correction for left-hand side of cell
        Mmn(q,b) = state(2 : numKpts, q)' * state(1 : numKpts - 1, neighbours.Nearest(q,b));
        
    case '2D'
        N = sqrt(numQpts); ind = 1;
        % Functions to translate the sides/corners by a reciprocal lattice vector
        [ds] = Sides2d(state, meshInd, recip.Size, recip.GkiIn);
        [dc] = Corners2d(state, meshInd, recip.Size, recip.GkiIn);
        % Corrections for neighbours in the reciprocal lattice basis directions
        b = 1; % 1st nearest neighbour vector, G1
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = meshInd(:,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,:,s), 1));
            ind = ind + 1;
        end
        b = 2; % 2nd nearest neighbour vector, -G1
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = meshInd(:,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,:,s),1));
            ind = ind + 1;
        end
        b = 3; % 3rd nearest neighbour vector, G2
        if sum(neighbours.WeightInd == b) > 0
            s = 3; % Correction for top of cell, G2
            q = meshInd(N,:);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,:,s), 1));
            ind = ind + 1;
        end
        b = 4; % 4th nearest neighbour vector, -G2
        if sum(neighbours.WeightInd == b) > 0
            s = 4; % Correction for bottom of cell, -G2
            q = meshInd(1,:);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,:,s), 1));
            ind = ind + 1;
        end
        % Corrections for diagonal neighbours
        b = 5; % 5th nearest neighbour vector, G1+G2
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = meshInd(1:N-1,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,2:N,s), 1));
            s = 3; % Correction for top of cell, G2
            q = meshInd(N,1:N-1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,2:N,s), 1));
            ci = 1; % Correction for top right-hand corner of cell, G1+G2
            q = meshInd(N,N);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
            ind = ind + 1;
        end
        b = 6; % 6th nearest neighbour vector, -(G1+G2)
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = meshInd(2:N,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,1:N-1,s), 1));
            s = 4; % Correction for bottom of cell, -G2
            q = meshInd(1,2:N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,1:N-1,s), 1));
            ci = 2; % Correction for bottom left-hand corner of cell, -(G1+G2)
            q = meshInd(1,1);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
            ind = ind + 1;
        end
        b = 7; % 7th nearest neighbour vector, G1-G2
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = meshInd(2:N,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,1:N-1,s), 1));
            s = 4; % Correction for bottom of cell, -G2
            q = meshInd(1,1:N-1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,2:N,s), 1));
            ci = 3; % Correction for bottom right-hand corner of cell, G1-G2
            q = meshInd(1,N);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
            ind = ind + 1;
        end
        b = 8; % 8th nearest neighbour vector, -(G1-G2)
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = meshInd(1:N-1,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,2:N,s), 1));
            s = 3; % Correction for top of cell, G2
            q = meshInd(N,2:N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* ds(:,1:N-1,s), 1));
            ci = 4; % Correction for top left-hand corner of cell, -(G1-G2)
            q = meshInd(N,1);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
        end
        
    case '3D'
        N = round(numQpts ^ (1/3)); ind = 1;
        % Functions to translate the sides/edegs/corners by a reciprocal lattice vector
        [ds] = Sides3d(state, meshInd, recip.Size, recip.GkiIn);
        [de] = Edges3d(state, meshInd, recip.Size, recip.GkiIn);
        [dc] = Corners3d(state, meshInd, recip.Size, recip.GkiIn);
        % Corrections for neighbours in the reciprocal lattice basis directions
        b = 1; % 1st nearest neighbour vector, G1
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = reshape(meshInd(:,N,:), 1, N ^ 2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,:,s), numKpts, N ^ 2), 1));
            ind = ind + 1;
        end
        b = 2; % 2nd nearest neighbour vector, -G1
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = reshape(meshInd(:,1,:), 1, N ^ 2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,:,s), numKpts, N ^ 2), 1));
            ind = ind + 1;
        end
        b = 3; % 3rd nearest neighbour vector, G2
        if sum(neighbours.WeightInd == b) > 0
            s = 3; % Correction for back of cell, G2
            q = reshape(meshInd(N,:,:), 1, N ^ 2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,:,s), numKpts, N ^ 2), 1));
            ind = ind + 1;
        end
        b = 4; % 4th nearest neighbour vector, -G2
        if sum(neighbours.WeightInd == b) > 0
            s = 4; % Correction for front of cell, -G2
            q = reshape(meshInd(1,:,:), 1, N ^ 2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,:,s), numKpts, N ^ 2), 1));
            ind = ind + 1;
        end
        b = 5; % 5th nearest neighbour vector, G3
        if sum(neighbours.WeightInd == b) > 0
            s = 5; % Correction for top of cell, G3
            q = reshape(meshInd(:,:,N), 1, N ^ 2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,:,s), numKpts, N ^ 2), 1));
            ind = ind + 1;
        end
        b = 6; % 6th nearest neighbour vector, -G3
        if sum(neighbours.WeightInd == b) > 0
            s = 6; % Correction for bottom of cell, -G3
            q = reshape(meshInd(:,:,1), 1, N ^ 2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,:,s), numKpts, N ^ 2), 1));
            ind = ind + 1;
        end
        % Corrections for diagonal neighbours in plane G1,G2
        b = 7; % 7th nearest neighbour vector, G1+G2
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = reshape(meshInd(1:N-1,N,:), 1, N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,:,s), numKpts, N*(N-1)),1));
            s = 3; % Correction for back of cell, G2
            q = reshape(meshInd(N,1:N-1,:),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,:,s), numKpts, N*(N-1)),1));
            e = 1; % Correction for back-right edge of cell, G1+G2
            q = meshInd(N,N,:);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        b = 8; % 8th nearest neighbour vector, -(G1+G2)
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = reshape(meshInd(2:N,1,:),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,:,s), numKpts, N*(N-1)),1));
            s = 4; % Correction for front of cell, -G2
            q = reshape(meshInd(1,2:N,:),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,:,s), numKpts, N*(N-1)),1));
            e = 2; % Correction for front-left edge of cell, -(G1+G2)
            q = meshInd(1,1,:);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        b = 9; % 9th nearest neighbour vector, G1-G2
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = reshape(meshInd(2:N,N,:),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,:,s), numKpts, N*(N-1)),1));
            s = 4; % Correction for front of cell, -G2
            q = reshape(meshInd(1,1:N-1,:),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,:,s), numKpts, N*(N-1)),1));
            e = 3; % Correction for front-right edge of cell, G1-G2
            q = meshInd(1,N,:);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        b = 10; % 10th nearest neighbour vector, -(G1-G2)
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = reshape(meshInd(1:N-1,1,:),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,:,s), numKpts, N*(N-1)),1));
            s = 3; % Correction for back of cell, G2
            q = reshape(meshInd(N,2:N,:),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,:,s), numKpts, N*(N-1)),1));
            e = 4; % Correction for back-left edge of cell, -(G1-G2)
            q = meshInd(N,1,:);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        % Corrections for diagonal neighbours in plane G1,G3
        b = 11; % 11th nearest neighbour vector, G1+G3
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = reshape(meshInd(:,N,1:N-1),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,2:N,s), numKpts, N*(N-1)),1));
            s = 5; % Correction for top of cell, G3
            q = reshape(meshInd(:,1:N-1,N),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,2:N,s), numKpts, N*(N-1)),1));
            e = 5; % Correction for top-right edge of cell, G1+G3
            q = meshInd(:,N,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        b = 12; % 12th nearest neighbour vector, -(G1+G3)
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = reshape(meshInd(:,1,2:N),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,1:N-1,s), numKpts, N*(N-1)),1));
            s = 6; % Correction for bottom of cell, -G3
            q = reshape(meshInd(:,2:N,1),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,1:N-1,s), numKpts, N*(N-1)),1));
            e = 6; % Correction for bottom-left edge of cell, -(G1+G3)
            q = meshInd(:,1,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        b = 13; % 13th nearest neighbour vector, G1-G3
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = reshape(meshInd(:,N,2:N),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,1:N-1,s), numKpts, N*(N-1)),1));
            s = 6; % Correction for bottom of cell, -G3
            q = reshape(meshInd(:,1:N-1,1),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,2:N,s), numKpts, N*(N-1)),1));
            e = 7; % Correction for bottom-right edge of cell, G1-G3
            q = meshInd(:,N,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        b = 14; % 14th nearest neighbour vector, -(G1-G3)
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = reshape(meshInd(:,1,1:N-1),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,2:N,s), numKpts, N*(N-1)),1));
            s = 5; % Correction for top of cell, G3
            q = reshape(meshInd(:,2:N,N),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,1:N-1,s), numKpts, N*(N-1)),1));
            e = 8; % Correction for top-left edge of cell, -(G1-G3)
            q = meshInd(:,1,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        % Corrections for diagonal neighbours in plane G2,G3
        b = 15; % 15th nearest neighbour vector, G2+G3
        if sum(neighbours.WeightInd == b) > 0
            s = 3; % Correction for back of cell, G2
            q = reshape(meshInd(N,:,1:N-1),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,2:N,s), numKpts, N*(N-1)),1));
            s = 5; % Correction for top of cell, G3
            q = reshape(meshInd(1:N-1,:,N),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,:,s), numKpts, N*(N-1)),1));
            e = 9; % Correction for top-back edge of cell, G2+G3
            q = meshInd(N,:,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        b = 16; % 16th nearest neighbour vector, -(G2+G3)
        if sum(neighbours.WeightInd == b) > 0
            s = 4; % Correction for front of cell, -G2
            q = reshape(meshInd(1,:,2:N),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,1:N-1,s), numKpts, N*(N-1)),1));
            s = 6; % Correction for bottom of cell, -G3
            q = reshape(meshInd(2:N,:,1),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,:,s), numKpts, N*(N-1)),1));
            e = 10; % Correction for bottom-front of cell, -(G2+G3)
            q = meshInd(1,:,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        b = 17; % 17th nearest neighbour vector, G2-G3
        if sum(neighbours.WeightInd == b) > 0
            s = 3; % Correction for back of cell, G2
            q = reshape(meshInd(N,:,2:N),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,1:N-1,s), numKpts, N*(N-1)),1));
            s = 6; % Correction for bottom of cell, -G3
            q = reshape(meshInd(1:N-1,:,1),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,:,s), numKpts, N*(N-1)),1));
            e = 11; % Correction for bottom-back edge of cell, G2-G3
            q = meshInd(N,:,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        b = 18; % 18th nearest neighbour vector, -(G2-G3)
        if sum(neighbours.WeightInd == b) > 0
            s = 4; % Correction for front of cell, -G2
            q = reshape(meshInd(1,:,1:N-1),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,:,2:N,s), numKpts, N*(N-1)),1));
            s = 5; % Correction for top of cell, G3
            q = reshape(meshInd(2:N,:,N),1,N*(N-1));
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,:,s), numKpts, N*(N-1)),1));
            e = 12; % Correction for top-front edge of cell, -(G2-G3)
            q = meshInd(1,:,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,:,e), 1));
            ind = ind + 1;
        end
        % Corrections for diagonal neighbours in plane G1+G2,G3
        b = 19; % 19th nearest neighbour vector, G1+G2+G3
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = reshape(meshInd(1:N-1,N,1:N-1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,2:N,s), numKpts, (N-1)^2), 1));
            s = 3; % Correction for back of cell, G2
            q = reshape(meshInd(N,1:N-1,1:N-1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,2:N,s), numKpts, (N-1)^2), 1));
            s = 5; % Correction for top of cell, G3
            q = reshape(meshInd(1:N-1,1:N-1,N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,2:N,s), numKpts, (N-1)^2), 1));
            e = 1; % Correction for back-right edge of cell, G1+G2
            q = meshInd(N,N,1:N-1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            e = 5; % Correction for top-right edge of cell, G1+G3
            q = meshInd(1:N-1,N,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            e = 9; % Correction for top-back edge of cell, G2+G3
            q = meshInd(N,1:N-1,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            ci = 1; % Correction for top-back-right corner of cell, G1+G2+G3
            q = meshInd(N,N,N);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
            ind = ind + 1;
        end
        b = 20; % 20th nearest neighbour vector, -(G1+G2+G3)
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, G1
            q = reshape(meshInd(2:N,1,2:N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,1:N-1,s), numKpts, (N-1)^2), 1));
            s = 4; % Correction for front of cell, G2
            q = reshape(meshInd(1,2:N,2:N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,1:N-1,s), numKpts, (N-1)^2), 1));
            s = 6; % Correction for bottom of cell, G3
            q = reshape(meshInd(2:N,2:N,1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,1:N-1,s), numKpts, (N-1)^2), 1));
            e = 2; % Correction for front-left edge of cell, -(G1+G2)
            q = meshInd(1,1,2:N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            e = 6; % Correction for bottom-left edge of cell, -(G1+G3)
            q = meshInd(2:N,1,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            e = 10; % Correction for bottom-front edge of cell, -(G2+G3)
            q = meshInd(1,2:N,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            ci = 2; % Correction for bottom-front-left corner of cell, -(G1+G2+G3)
            q = meshInd(1,1,1);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
            ind = ind + 1;
        end
        b = 21; % 21st nearest neighbour vector, G1+G2-G3
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = reshape(meshInd(1:N-1,N,2:N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,1:N-1,s), numKpts, (N-1)^2), 1));
            s = 3; % Correction for back of cell, G2
            q = reshape(meshInd(N,1:N-1,2:N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,1:N-1,s), numKpts, (N-1)^2), 1));
            s = 6; % Correction for bottom of cell, -G3
            q = reshape(meshInd(1:N-1,1:N-1,1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)).*reshape(ds(:,2:N,2:N,s), numKpts, (N-1)^2), 1));
            e = 1; % Correction for back-right edge of cell, G1+G2
            q = meshInd(N,N,2:N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            e = 7; % Correction for bottom-right edge of cell, G1-G3
            q = meshInd(1:N-1,N,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            e = 11; % Correction for bottom-back edge of cell, G2-G3
            q = meshInd(N,1:N-1,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            ci = 3; % Correction for bottom-back-right corner of cell, G1+G2-G3
            q = meshInd(N,N,1);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
            ind = ind + 1;
        end
        b = 22; % 22nd nearest neighbour vector, -(G1+G2-G3)
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = reshape(meshInd(2:N,1,1:N-1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,2:N,s), numKpts, (N-1)^2), 1));
            s = 4; % Correction for front of cell, -G2
            q = reshape(meshInd(1,2:N,1:N-1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,2:N,s), numKpts, (N-1)^2), 1));
            s = 5; % Correction for top of cell, G3
            q = reshape(meshInd(2:N,2:N,N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,1:N-1,s), numKpts, (N-1)^2), 1));
            e = 2; % Correction for front-left edge of cell, -(G1+G2)
            q = meshInd(1,1,1:N-1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            e = 8; % Correction for top-left edge of cell, -(G1-G3)
            q = meshInd(2:N,1,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            e = 12; % Correction for top-front edge of cell, -(G2-G3)
            q = meshInd(1,2:N,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            ci = 4; % Correction for top-front-left corner of cell, -(G1+G2-G3)
            q = meshInd(1,1,N);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
            ind = ind + 1;
        end
        % Corrections for diagonal neighbours in plane G1-G2,G3
        b = 23; % 23rd nearest neighbour vector, G1-G2+G3
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = reshape(meshInd(2:N,N,1:N-1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,2:N,s), numKpts, (N-1)^2), 1));
            s = 4; % Correction for front of cell, -G2
            q = reshape(meshInd(1,1:N-1,1:N-1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,2:N,s), numKpts, (N-1)^2), 1));
            s = 5; % Correction for top of cell, G3
            q = reshape(meshInd(2:N,1:N-1,N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,2:N,s), numKpts, (N-1)^2), 1));
            e = 3; % Correction for front-right edge of cell, G1-G2
            q = meshInd(1,N,1:N-1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            e = 5; % Correction for top-right edge of cell, G1+G3
            q = meshInd(2:N,N,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            e = 12; % Correction for top-front edge of cell, -G2+G3
            q = meshInd(1,1:N-1,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            ci = 5; % Correction for top-front-right corner of cell, G1-G2+G3
            q = meshInd(1,N,N);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
            ind = ind + 1;
        end
        b = 24; % 24th nearest neighbour vector, -(G1-G2+G3)
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = reshape(meshInd(1:N-1,1,2:N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,1:N-1,s), numKpts, (N-1)^2), 1));
            s = 3; % Correction for back of cell, G2
            q = reshape(meshInd(N,2:N,2:N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,1:N-1,s), numKpts, (N-1)^2), 1));
            s = 6; % Correction for bottom of cell, -G3
            q = reshape(meshInd(1:N-1,2:N,1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,1:N-1,s), numKpts, (N-1)^2), 1));
            e = 4; % Correction for back-left edge of cell, -(G1-G2)
            q = meshInd(N,1,2:N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            e = 6; % Correction for bottom-left edge of cell, -(G1+G3)
            q = meshInd(1:N-1,1,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            e = 11; % Correction for bottom-back edge of cell, G2-G3
            q = meshInd(N,2:N,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            ci = 6; % Correction for bottom-back-left corner of cell, -(G1-G2+G3)
            q = meshInd(N,1,1);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
            ind = ind + 1;
        end
        b = 25; % 25th nearest neighbour vector, G1-G2-G3
        if sum(neighbours.WeightInd == b) > 0
            s = 1; % Correction for right-hand side of cell, G1
            q = reshape(meshInd(2:N,N,2:N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,1:N-1,s), numKpts, (N-1)^2), 1));
            s = 4; % Correction for front of cell, -G2
            q = reshape(meshInd(1,1:N-1,2:N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,1:N-1,s), numKpts, (N-1)^2), 1));
            s = 6; % Correction for bottom of cell, -G3
            q = reshape(meshInd(2:N,1:N-1,1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,2:N,s), numKpts, (N-1)^2), 1));
            e = 3; % Correction for front-right edge of cell, G1-G2
            q = meshInd(1,N,2:N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            e = 7; % Correction for bottom-right edge of cell, G1-G3
            q = meshInd(2:N,N,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            e = 10; % Correction for bottom-front edge of cell, -(G2+G3)
            q = meshInd(1,1:N-1,1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            ci = 7; % Correction for bottom-front-right corner of cell, G1-G2-G3
            q = meshInd(1,N,1);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
            ind = ind + 1;
        end
        b = 26; % 26th nearest neighbour vector, -(G1-G2-G3)
        if sum(neighbours.WeightInd == b) > 0
            s = 2; % Correction for left-hand side of cell, -G1
            q = reshape(meshInd(1:N-1,1,1:N-1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,2:N,s), numKpts, (N-1)^2), 1));
            s = 3; % Correction for back of cell, G2
            q = reshape(meshInd(N,2:N,1:N-1),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,1:N-1,2:N,s), numKpts, (N-1)^2), 1));
            s = 5; % Correction for top of cell, G3
            q = reshape(meshInd(1:N-1,2:N,N),1,(N-1)^2);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* reshape(ds(:,2:N,1:N-1,s), numKpts, (N-1)^2), 1));
            e = 4; % Correction for back-left edge of cell, -G1+G2
            q = meshInd(N,1,1:N-1);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            e = 8; % Correction for top-left edge of cell, -G1+G3
            q = meshInd(1:N-1,1,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,2:N,e), 1));
            e = 9; % Correction for top-back edge of cell, G2+G3
            q = meshInd(N,2:N,N);
            Mmn(q, ind) = transpose(sum(conj(state(:,q)) .* de(:,1:N-1,e), 1));
            ci = 8; % Correction for top-back-left corner of cell, -G1+G2+G3
            q = meshInd(N,1,N);
            Mmn(q, ind) = state(:,q)' * dc(:,ci);
        end
end

% -----------------------------------------------------------------------------------------------------------
% Translates the Bloch states on each side of a 2D cell by one reciprocal lattice vector to the other side of
% the cell. s labels the sides.
% -----------------------------------------------------------------------------------------------------------
function  [ds] = Sides2d(c, meshInd, D, gkiIn)
numQpts = size(c,2); N = sqrt(numQpts);
cTemp = zeros(prod(D), numQpts); cTemp(gkiIn, :) = c;
cTemp = reshape(cTemp, D(2), D(1), numQpts);
ds = zeros(D(2), D(1), N, 4);
% s = 1; G1; Shift the states on the left-hand side of the cell one reciprocal lattice vector to the right
q = meshInd(:,1); ds(:,1:D(1)-1,:,1) = cTemp(:,2:D(1),q);
% s = 2; -G1; Shift the states on the right-hand side of the cell one reciprocal lattice vector to the left
q = meshInd(:,N); ds(:,2:D(1),:,2) = cTemp(:,1:D(1)-1,q);
% s = 3; G2; Shift the states on the bottom of the cell one reciprocal lattice vector upwards
q = meshInd(1,:); ds(1:D(2)-1,:,:,3) = cTemp(2:D(2),:,q);
% s = 4; -G2; Shift the states on the top of the cell one reciprocal lattice vector downwards
q = meshInd(N,:); ds(2:D(2),:,:,4) = cTemp(1:D(2)-1,:,q);
% Reshape ds
ds = reshape(ds, prod(D), N, 4); ds = ds(gkiIn,:,:);

% --------------------------------------------------------------------------------------------------------
% Translates the Bloch states on each corner of a 2D cell by one reciprocal lattice vector to the opposite
% corner of the cell. ci labels the corners.
% --------------------------------------------------------------------------------------------------------
function  [dc] = Corners2d(c, meshInd, D, gkiIn)
numQpts = size(c,2); N = sqrt(numQpts);
cTemp = zeros(prod(D), numQpts); cTemp(gkiIn,:) = c;
cTemp = reshape(cTemp, D(2), D(1), numQpts);
dc = zeros(D(2), D(1), 4);
% ci = 1; G1+G2; Shift the states on the bottom left-hand corner of the cell one reciprocal lattice vector to
% the top-right
dc(1:D(2)-1,1:D(1)-1,1) = cTemp(2:D(2),2:D(1),meshInd(1,1));
% ci = 2; -(G1+G2); Shift the states on the top right-hand corner of the cell one reciprocal lattice vector to
% the bottom-left
dc(2:D(2),2:D(1),2) = cTemp(1:D(2)-1,1:D(1)-1,meshInd(N,N));
% ci = 3; G1-G2; Shift the states on the top left-hand corner of the cell one reciprocal lattice vector to the
% bottom-right
dc(2:D(2),1:D(1)-1,3) = cTemp(1:D(2)-1,2:D(1),meshInd(N,1));
% ci = 4; -(G1-G2); Shift the states on the bottom right-hand corner of the cell one reciprocal lattice vector
% to the top-left
dc(1:D(2)-1,2:D(1),4) = cTemp(2:D(2),1:D(1)-1,meshInd(1,N));
% Reshape dc
dc = reshape(dc, prod(D), 4); dc = dc(gkiIn,:);

% -----------------------------------------------------------------------------------------------------------
% Translates the Bloch states on each side of a 3D cell by one reciprocal lattice vector to the other side of
% the cell. s labels the sides.
% -----------------------------------------------------------------------------------------------------------
function  [ds] = Sides3d(c, meshInd, D, gkiIn)
numQpts = size(c,2); N = round(numQpts ^ (1/3));
cTemp = zeros(prod(D), numQpts); cTemp(gkiIn,:) = c;
cTemp = reshape(cTemp, D(2), D(1), D(3), numQpts);
ds = zeros(D(2), D(1), D(3), N ^ 2, 6);
% s = 1; G1; Shift the states on the left-hand side of the cell one reciprocal lattice vector to the right
q = reshape(meshInd(:,1,:),1,N^2); ds(:,1:D(1)-1,:,:,1) = cTemp(:,2:D(1),:,q);
% s = 2; -G1; Shift the states on the right-hand side of the cell one reciprocal lattice vector to the left
q = reshape(meshInd(:,N,:),1,N^2); ds(:,2:D(1),:,:,2) = cTemp(:,1:D(1)-1,:,q);
% s = 3; G2; Shift the states on the front of the cell one reciprocal lattice vector backwards
q = reshape(meshInd(1,:,:),1,N^2); ds(1:D(2)-1,:,:,:,3) = cTemp(2:D(2),:,:,q);
% s = 4; -G2; Shift the states on the back of the cell one reciprocal lattice vector forwards
q = reshape(meshInd(N,:,:),1,N^2); ds(2:D(2),:,:,:,4) = cTemp(1:D(2)-1,:,:,q);
% s = 5; G3; Shift the states on the bottom of the cell one reciprocal lattice vector upwards
q = reshape(meshInd(:,:,1),1,N^2); ds(:,:,1:D(3)-1,:,5) = cTemp(:,:,2:D(3),q);
% s = 6; -G3; Shift the states on the top of the cell one reciprocal lattice vector downwards
q = reshape(meshInd(:,:,N),1,N^2); ds(:,:,2:D(3),:,6) = cTemp(:,:,1:D(3)-1,q);
% Reshape ds
ds = reshape(ds, prod(D), N, N, 6); ds = ds(gkiIn,:,:,:);

% -----------------------------------------------------------------------------------------------------------
% Translates the Bloch states on each edge of a 3D cell by one reciprocal lattice vector to the opposite edge
% of the cell. e labels the edges.
% -----------------------------------------------------------------------------------------------------------
function  [de] = Edges3d(c, meshInd, D, gkiIn)
numQpts = size(c,2); N = round(numQpts ^ (1/3));
cTemp = zeros(prod(D), numQpts); cTemp(gkiIn,:) = c;
cTemp = reshape(cTemp, D(2), D(1), D(3), numQpts);
de = zeros(D(2), D(1), D(3), N, 12);
% e = 1; G1+G2; Shift the states on the front-left edge of the cell one reciprocal lattice vector to the
% back-right
de(1:D(2)-1,1:D(1)-1,:,:,1) = cTemp(2:D(2),2:D(1),:,meshInd(1,1,:));
% e = 2; -(G1+G2); Shift the states on the back-right edge of the cell one reciprocal lattice vector to the
% front-left
de(2:D(2),2:D(1),:,:,2) = cTemp(1:D(2)-1,1:D(1)-1,:,meshInd(N,N,:));
% e = 3; G1-G2; Shift the states on the back-left edge of the cell one reciprocal lattice vector to the
% front-right
de(2:D(2),1:D(1)-1,:,:,3) = cTemp(1:D(2)-1,2:D(1),:,meshInd(N,1,:));
% e = 4; -(G1-G2); Shift the states on the front-right edge of the cell one reciprocal lattice vector to the
% back-left
de(1:D(2)-1,2:D(1),:,:,4) = cTemp(2:D(2),1:D(1)-1,:,meshInd(1,N,:));
% e = 5; G1+G3; Shift the states on the bottom-left edge of the cell one reciprocal lattice vector to the
% top-right
de(:,1:D(1)-1,1:D(3)-1,:,5) = cTemp(:,2:D(1),2:D(3),meshInd(:,1,1));
% e = 6; -(G1+G3); Shift the states on the top-right edge of the cell one reciprocal lattice vector to the
% bottom-left
de(:,2:D(1),2:D(3),:,6) = cTemp(:,1:D(1)-1,1:D(3)-1,meshInd(:,N,N));
% e = 7; G1-G3; Shift the states on the top-left edge of the cell one reciprocal lattice vector to the
% bottom-right
de(:,1:D(1)-1,2:D(3),:,7) = cTemp(:,2:D(1),1:D(3)-1,meshInd(:,1,N));
% e = 8; -(G1-G3); Shift the states on the bottom-right edge of the cell one reciprocal lattice vector to the
% top-left
de(:,2:D(1),1:D(3)-1,:,8) = cTemp(:,1:D(1)-1,2:D(3),meshInd(:,N,1));
% e = 9; G2+G3; Shift the states on the bottom-front edge of the cell one reciprocal lattice vector to the
% top-back
de(1:D(2)-1,:,1:D(3)-1,:,9) = cTemp(2:D(2),:,2:D(3),meshInd(1,:,1));
% e = 10; -(G2+G3); Shift the states on the top-back edge of the cell one reciprocal lattice vector to the
% bottom-front
de(2:D(2),:,2:D(3),:,10) = cTemp(1:D(2)-1,:,1:D(3)-1,meshInd(N,:,N));
% e = 11; G2-G3; Shift the states on the top-front edge of the cell one reciprocal lattice vector to the
% bottom-back
de(1:D(2)-1,:,2:D(3),:,11) = cTemp(2:D(2),:,1:D(3)-1,meshInd(1,:,N));
% e = 12; -(G2-G3); Shift the states on the bottom-back edge of the cell one reciprocal lattice vector to the
% top-front
de(2:D(2),:,1:D(3)-1,:,12) = cTemp(1:D(2)-1,:,2:D(3),meshInd(N,:,1));
% Reshape de
de = reshape(de, prod(D), N, 12); de = de(gkiIn,:,:);

% --------------------------------------------------------------------------------------------------------
% Translates the Bloch states on each corner of a 3D cell by one reciprocal lattice vector to the opposite
% corner of the cell. ci labels the corners.
% --------------------------------------------------------------------------------------------------------
function  [dc] = Corners3d(c, meshInd, D, gkiIn)
numQpts = size(c,2); N = round(numQpts ^ (1/3));
cTemp = zeros(prod(D), numQpts); cTemp(gkiIn,:) = c;
cTemp = reshape(cTemp, D(2), D(1), D(3), numQpts);
dc = zeros(D(2), D(1), D(3), 8);
% ci = 1; G1+G2+G3; Shift the states on the bottom-front-left corner of the cell one reciprocal lattice vector
% to the top-back-right
dc(1:D(2)-1,1:D(1)-1,1:D(3)-1,1) = cTemp(2:D(2),2:D(1),2:D(3),meshInd(1,1,1));
% ci = 2; -(G1+G2+G3); Shift the states on the top-back-right corner of the cell one reciprocal lattice vector
% to the bottom-front-left
dc(2:D(2),2:D(1),2:D(3),2) = cTemp(1:D(2)-1,1:D(1)-1,1:D(3)-1,meshInd(N,N,N));
% ci = 3; G1+G2-G3; Shift the states on the top-front-left corner of the cell one reciprocal lattice vector to
% the bottom-back-right
dc(1:D(2)-1,1:D(1)-1,2:D(3),3) = cTemp(2:D(2),2:D(1),1:D(3)-1,meshInd(1,1,N));
% ci = 4; -(G1+G2-G3); Shift the states on the bottom-back-right corner of the cell one reciprocal lattice
% vector to the top-front-left
dc(2:D(2),2:D(1),1:D(3)-1,4) = cTemp(1:D(2)-1,1:D(1)-1,2:D(3),meshInd(N,N,1));
% ci = 5; G1-G2+G3; Shift the states on the bottom-back-left corner of the cell one reciprocal lattice vector
% to the top-front-right
dc(2:D(2),1:D(1)-1,1:D(3)-1,5) = cTemp(1:D(2)-1,2:D(1),2:D(3),meshInd(N,1,1));
% ci = 6; -(G1-G2+G3); Shift the states on the top-front-right corner of the cell one reciprocal lattice
% vector to the bottom-back-left
dc(1:D(2)-1,2:D(1),2:D(3),6) = cTemp(2:D(2),1:D(1)-1,1:D(3)-1,meshInd(1,N,N));
% ci = 7; G1-G2-G3; Shift the states on the top-back-left corner of the cell one reciprocal lattice vector to
% the bottom-front-right
dc(2:D(2),1:D(1)-1,2:D(3),7) = cTemp(1:D(2)-1,2:D(1),1:D(3)-1,meshInd(N,1,N));
% ci = 8; -(G1-G2-G3); Shift the states on the bottom-front-right corner of the cell one reciprocal lattice
% vector to the top-back-left
dc(1:D(2)-1,2:D(1),1:D(3)-1,8) = cTemp(2:D(2),1:D(1)-1,2:D(3),meshInd(1,N,1));
% Reshape dc
dc = reshape(dc, prod(D), 8); dc = dc(gkiIn,:);
