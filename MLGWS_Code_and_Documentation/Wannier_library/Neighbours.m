
% Determine the nearest neighbours and associated weights for a Bloch mesh
%
% PROPERTIES
% B: Vectors to each nearest neighbour mesh point
% Weight: Weight corresponding to each B vector
% WeightInd: Index of the B vector to each neighbour
% Nearest: Mesh index of the nearest neighbour, labelled by q,b
%
% METHODS
% Neighbours(lattice, N, meshInd): Constructor. 'mesh' must be a bloch.QMesh instance of SuperCell.
%
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

classdef Neighbours
    
    properties(GetAccess = 'public', SetAccess = 'private')
        % DECLARE PROPERTIES AND DEFAULT VALUES
        B;
        Weight;
        WeightInd;
        Nearest;
    end
    
    methods
        function neighbours = Neighbours(lattice, N, meshInd)
            % CLASS CONSTRUCTOR
            %
            % REQUIRED INPUT ARGUMENTS
            % lattice --> instance of the Lattice class
            %       N --> number of mesh points in each direction
            % meshInd --> index of the mesh points
            
            switch lattice.Dimension
                case '1D'
                    b(:,1) = lattice.G(:,1);
                    b(:,2) = -b(:,1);
                case '2D'
                    b(:,1) = lattice.G(:,1);
                    b(:,2) = -b(:,1);
                    b(:,3) = lattice.G(:,2);
                    b(:,4) = -b(:,3);
                    b(:,5) = lattice.G(:,1) + lattice.G(:,2);
                    b(:,6) = -b(:,5);
                    b(:,7) = lattice.G(:,1) - lattice.G(:,2);
                    b(:,8) = -b(:,7);
                case '3D'
                    b(:,1) = lattice.G(:,1);
                    b(:,2) = -b(:,1);
                    b(:,3) = lattice.G(:,2);
                    b(:,4) = -b(:,3);
                    b(:,5) = lattice.G(:,3);
                    b(:,6) = -b(:,5);
                    b(:,7) = lattice.G(:,1) + lattice.G(:,2);
                    b(:,8) = -b(:,7);
                    b(:,9) = lattice.G(:,1) - lattice.G(:,2);
                    b(:,10) = -b(:,9);
                    b(:,11) = lattice.G(:,1) + lattice.G(:,3);
                    b(:,12) = -b(:,11);
                    b(:,13) = lattice.G(:,1) - lattice.G(:,3);
                    b(:,14) = -b(:,13);
                    b(:,15) = lattice.G(:,2) + lattice.G(:,3);
                    b(:,16) = -b(:,15);
                    b(:,17) = lattice.G(:,2) - lattice.G(:,3);
                    b(:,18) = -b(:,17);
                    b(:,19) = lattice.G(:,1) + lattice.G(:,2) + lattice.G(:,3);
                    b(:,20) = -b(:,19);
                    b(:,21) = lattice.G(:,1) + lattice.G(:,2) - lattice.G(:,3);
                    b(:,22) = -b(:,21);
                    b(:,23) = lattice.G(:,1) - lattice.G(:,2) + lattice.G(:,3);
                    b(:,24) = -b(:,23);
                    b(:,25) = lattice.G(:,1) - lattice.G(:,2) - lattice.G(:,3);
                    b(:,26) = -b(:,25);
            end
            b = b / N;
            
% ------------------------------------------------------------------------------------------------------------
% Calculate the weights associated with each neighbour. The weights must satisfy A*w = qw, where qw = [1;1;0].
% See Sec. 3.2 of Ref. 2. Iterate through shells of increasing length until the condition is satisfied.
% ------------------------------------------------------------------------------------------------------------
            switch lattice.Dimension
                case '1D'
                    wb(1:2) = 1 / (2 * b(1,1) ^ 2); % Trivial calculation in 1D
                    wbInd = 0; % Unused variable in 1D
                case {'2D', '3D'}
                    bNorm = sqrt(sum(b .^ 2, 1)); [bNorm, bInd] = sort(bNorm); % Lengths of the B-vectors
                    % Initial values
                    if strcmp(lattice.Dimension, '2D')
                        qw = [1; 1; 0]; qwTest = [0; 0; 0];
                    elseif strcmp(lattice.Dimension, '3D')
                        qw = [1; 1; 1; 0; 0; 0]; qwTest = [0; 0; 0; 0; 0; 0];
                    end
                    loop = 0; s = 1; wbInd = [];
                    % Iterate through shells of increasing length until qw_test == qw
                    while max(abs(qwTest - qw)) > 1e-9
                        wbIndS = bInd(bNorm == bNorm(loop + 1)); % Index of the B vectors in shell s
                        numB(s) = length(wbIndS); % No. of B vectors in shell s
                        wbInd = [wbInd wbIndS]; % Index of all B vectors, sorted by weight, up to and including shell s
                        bs = b(:, wbIndS); % B vectors in the shell s
                        % A(j,s) = sum over b of B_b(alpha) * B_b(beta) for shell s and j = alpha.beta
                        if strcmp(lattice.Dimension, '2D')
                            A(1,s) = sum(bs(1,:) .^ 2, 2); % xx
                            A(2,s) = sum(bs(2,:) .^ 2, 2); % yy
                            A(3,s) = sum(bs(1,:) .* bs(2,:,:), 2); % xy
                        elseif strcmp(lattice.Dimension, '3D')
                            A(1,s) = sum(bs(1,:,:) .^ 2, 2); % xx
                            A(2,s) = sum(bs(2,:,:) .^ 2, 2); % yy
                            A(3,s) = sum(bs(3,:,:) .^ 2, 2); % zz
                            A(4,s) = sum(bs(1,:,:) .* bs(2,:,:), 2); % xy
                            A(5,s) = sum(bs(1,:,:) .* bs(3,:,:), 2); % xz
                            A(6,s) = sum(bs(2,:,:) .* bs(3,:,:), 2); % yz
                        end
                        [U,S,V] = svd(A); % Calculate the SVD of A, Eq. 26 of Ref. 2
                        w = transpose(V / S * U' * qw); % Weight associcated to each shell s, Eq. 27 of Ref. 2
                        wb = [];
                        for s2 = 1 : s
                            wb = [wb w(s2 * ones(1, numB(s2)))]; end % Weight of each B vector
                        qwTest = A * transpose(w); % Check that A*w = q, Eq. 25 of Ref. 2
                        loop = loop + numB(s); s = s + 1;
                    end
                    [wbIndTemp, wbInd2] = sort(wbInd);
                    wb = wb(wbInd2); % Sort the weights in original order of B-vectors
                    wbInd = sort(wbInd); % Now the index of included B-vectors, sorted in original order
                    b = b(:, wbInd); % Remove neighbours with zero weight
            end
            neighbours.B = b;
            neighbours.Weight = wb;
            neighbours.WeightInd = wbInd;
            
% -------------------------------------------------
% Find the index of each neighbour for all q-points
% -------------------------------------------------
            switch lattice.Dimension
                case '1D'
                    % Find the index of each neighbour for all q-points
                    neighbours.Nearest = [2:N,1;N,1:N-1];
                    neighbours.Nearest = transpose(neighbours.Nearest); % q,b
                case '2D'
                    numBVecs = length(wb); numQpts = N ^ 2;
                    neighbours.Nearest = zeros(N,N,numBVecs); % y,x,b
                    % 1st B vector neighbours, G1
                    neighbours.Nearest(:,1:N-1,1) = meshInd(:,2:N); % Index body
                    neighbours.Nearest(:,N,1) = meshInd(:,1); % Index edge 1
                    % 2nd B vector neighbours, -G1
                    neighbours.Nearest(:,2:N,2) = meshInd(:,1:N-1); % Index body
                    neighbours.Nearest(:,1,2) = meshInd(:,N); % Index edge 1
                    % 3rd B vector neighbours, G2
                    neighbours.Nearest(1:N-1,:,3) = meshInd(2:N,:); % Index body
                    neighbours.Nearest(N,:,3) = meshInd(1,:); % Index edge 2
                    % 4th B vector neighbours, -G2
                    neighbours.Nearest(2:N,:,4) = meshInd(1:N-1,:); % Index body
                    neighbours.Nearest(1,:,4) = meshInd(N,:); % Index edge 2
                    % 5th B vector neighbours, G1+G2
                    neighbours.Nearest(1:N-1,1:N-1,5) = meshInd(2:N,2:N); % Index body
                    neighbours.Nearest(1:N-1,N,5) = meshInd(2:N,1); % Index edge 1
                    neighbours.Nearest(N,1:N-1,5) = meshInd(1,2:N); % Index edge 2
                    neighbours.Nearest(N,N,5) = meshInd(1,1); % Index corner
                    % 6th B vector neighbours, -(G1+G2)
                    neighbours.Nearest(2:N,2:N,6) = meshInd(1:N-1,1:N-1); % Index body
                    neighbours.Nearest(2:N,1,6) = meshInd(1:N-1,N); % Index edge 1
                    neighbours.Nearest(1,2:N,6) = meshInd(N,1:N-1); % Index edge 2
                    neighbours.Nearest(1,1,6) = meshInd(N,N); % Index corner
                    % 7th B vector neighbours, G1-G2
                    neighbours.Nearest(2:N,1:N-1,7) = meshInd(1:N-1,2:N); % Index body
                    neighbours.Nearest(2:N,N,7) = meshInd(1:N-1,1); % Index edge 1
                    neighbours.Nearest(1,1:N-1,7) = meshInd(N,2:N); % Index edge 2
                    neighbours.Nearest(1,N,7) = meshInd(N,1); % Index corner
                    % 8th B vector neighbours, -(G1-G2)
                    neighbours.Nearest(1:N-1,2:N,8) = meshInd(2:N,1:N-1); % Index body
                    neighbours.Nearest(1:N-1,1,8) = meshInd(2:N,N); % Index edge 1
                    neighbours.Nearest(N,2:N,8) = meshInd(1,1:N-1); % Index edge 2
                    neighbours.Nearest(N,1,8) = meshInd(1,N); % Index corner
                    % Remove neighbour indices with zero weight and reshape the nearest neighbour index
                    neighbours.Nearest = neighbours.Nearest(:,:,wbInd);
                    neighbours.Nearest = reshape(neighbours.Nearest,numQpts,numBVecs); % q,b
                case '3D'
                    numBVecs = length(wb); numQpts = N ^ 3;
                    neighbours.Nearest = zeros(N,N,N,12); % y,x,z,b
                    % 1st B vector neighbours, G1
                    neighbours.Nearest(:,1:N-1,:,1) = meshInd(:,2:N,:); % Index body
                    neighbours.Nearest(:,N,:,1) = meshInd(:,1,:); % Index face 1
                    % 2nd B vector neighbours, -G1
                    neighbours.Nearest(:,2:N,:,2) = meshInd(:,1:N-1,:); % Index body
                    neighbours.Nearest(:,1,:,2) = meshInd(:,N,:); % Index face 1
                    % 3rd B vector neighbours, G2
                    neighbours.Nearest(1:N-1,:,:,3) = meshInd(2:N,:,:); % Index body
                    neighbours.Nearest(N,:,:,3) = meshInd(1,:,:); % Index face 2
                    % 4th B vector neighbours, -G2
                    neighbours.Nearest(2:N,:,:,4) = meshInd(1:N-1,:,:); % Index body
                    neighbours.Nearest(1,:,:,4) = meshInd(N,:,:); % Index face 2
                    % 5th B vector neighbours, G3
                    neighbours.Nearest(:,:,1:N-1,5) = meshInd(:,:,2:N); % Index body
                    neighbours.Nearest(:,:,N,5) = meshInd(:,:,1); % Index face 3
                    % 6th B vector neighbours, -G3
                    neighbours.Nearest(:,:,2:N,6) = meshInd(:,:,1:N-1); % Index body
                    neighbours.Nearest(:,:,1,6) = meshInd(:,:,N); % Index face 3
                    % 7th B vector neighbours, G1+G2
                    neighbours.Nearest(1:N-1,1:N-1,:,7) = meshInd(2:N,2:N,:); % Index body
                    neighbours.Nearest(1:N-1,N,:,7) = meshInd(2:N,1,:); % Index face 1
                    neighbours.Nearest(N,1:N-1,:,7) = meshInd(1,2:N,:); % Index face 2
                    neighbours.Nearest(N,N,:,7) = meshInd(1,1,:); % Index edge 12
                    % 8th B vector neighbours, -(G1+G2)
                    neighbours.Nearest(2:N,2:N,:,8) = meshInd(1:N-1,1:N-1,:); % Index body
                    neighbours.Nearest(2:N,1,:,8) = meshInd(1:N-1,N,:); % Index face 1
                    neighbours.Nearest(1,2:N,:,8) = meshInd(N,1:N-1,:); % Index face 2
                    neighbours.Nearest(1,1,:,8) = meshInd(N,N,:); % Index edge 12
                    % 9th B vector neighbours, G1-G2
                    neighbours.Nearest(2:N,1:N-1,:,9) = meshInd(1:N-1,2:N,:); % Index body
                    neighbours.Nearest(2:N,N,:,9) = meshInd(1:N-1,1,:); % Index face 1
                    neighbours.Nearest(1,1:N-1,:,9) = meshInd(N,2:N,:); % Index face 2
                    neighbours.Nearest(1,N,:,9) = meshInd(N,1,:); % Index edge 12
                    % 10th B vector neighbours, -(G1-G2)
                    neighbours.Nearest(1:N-1,2:N,:,10) = meshInd(2:N,1:N-1,:); % Index body
                    neighbours.Nearest(1:N-1,1,:,10) = meshInd(2:N,N,:); % Index face 1
                    neighbours.Nearest(N,2:N,:,10) = meshInd(1,1:N-1,:); % Index face 2
                    neighbours.Nearest(N,1,:,10) = meshInd(1,N,:); % Index edge 12
                    % 11th B vector neighbours, G1+G3
                    neighbours.Nearest(:,1:N-1,1:N-1,11) = meshInd(:,2:N,2:N); % Index body
                    neighbours.Nearest(:,N,1:N-1,11) = meshInd(:,1,2:N); % Index face 1
                    neighbours.Nearest(:,1:N-1,N,11) = meshInd(:,2:N,1); % Index face 3
                    neighbours.Nearest(:,N,N,11) = meshInd(:,1,1); % Index edge 13
                    % 12th B vector neighbours, -(G1+G3)
                    neighbours.Nearest(:,2:N,2:N,12) = meshInd(:,1:N-1,1:N-1); % Index body
                    neighbours.Nearest(:,1,2:N,12) = meshInd(:,N,1:N-1); % Index face 1
                    neighbours.Nearest(:,2:N,1,12) = meshInd(:,1:N-1,N); % Index face 3
                    neighbours.Nearest(:,1,1,12) = meshInd(:,N,N); % Index edge 13
                    % 13th B vector neighbours, G1-G3
                    neighbours.Nearest(:,1:N-1,2:N,13) = meshInd(:,2:N,1:N-1); % Index body
                    neighbours.Nearest(:,N,2:N,13) = meshInd(:,1,1:N-1); % Index face 1
                    neighbours.Nearest(:,1:N-1,1,13) = meshInd(:,2:N,N); % Index face 3
                    neighbours.Nearest(:,N,1,13) = meshInd(:,1,N); % Index edge 13
                    % 14th B vector neighbours, -(G1-G3)
                    neighbours.Nearest(:,2:N,1:N-1,14) = meshInd(:,1:N-1,2:N); % Index body
                    neighbours.Nearest(:,1,1:N-1,14) = meshInd(:,N,2:N); % Index face 1
                    neighbours.Nearest(:,2:N,N,14) = meshInd(:,1:N-1,1); % Index face 3
                    neighbours.Nearest(:,1,N,14) = meshInd(:,N,1); % Index edge 13
                    % 15th B vector neighbours, G2+G3
                    neighbours.Nearest(1:N-1,:,1:N-1,15) = meshInd(2:N,:,2:N); % Index body
                    neighbours.Nearest(N,:,1:N-1,15) = meshInd(1,:,2:N); % Index face 2
                    neighbours.Nearest(1:N-1,:,N,15) = meshInd(2:N,:,1); % Index face 3
                    neighbours.Nearest(N,:,N,15) = meshInd(1,:,1); % Index edge 23
                    % 16th B vector neighbours, -(G2+G3)
                    neighbours.Nearest(2:N,:,2:N,16) = meshInd(1:N-1,:,1:N-1); % Index body
                    neighbours.Nearest(1,:,2:N,16) = meshInd(N,:,1:N-1); % Index face 2
                    neighbours.Nearest(2:N,:,1,16) = meshInd(1:N-1,:,N); % Index face 3
                    neighbours.Nearest(1,:,1,16) = meshInd(N,:,N); % Index edge 23
                    % 17th B vector neighbours, G2-G3
                    neighbours.Nearest(1:N-1,:,2:N,17) = meshInd(2:N,:,1:N-1); % Index body
                    neighbours.Nearest(N,:,2:N,17) = meshInd(1,:,1:N-1); % Index face 2
                    neighbours.Nearest(1:N-1,:,1,17) = meshInd(2:N,:,N); % Index face 3
                    neighbours.Nearest(N,:,1,17) = meshInd(1,:,N); % Index edge 23
                    % 18th B vector neighbours, -(G2-G3)
                    neighbours.Nearest(2:N,:,1:N-1,18) = meshInd(1:N-1,:,2:N); % Index body
                    neighbours.Nearest(1,:,1:N-1,18) = meshInd(N,:,2:N); % Index face 2
                    neighbours.Nearest(2:N,:,N,18) = meshInd(1:N-1,:,1); % Index face 3
                    neighbours.Nearest(1,:,N,18) = meshInd(N,:,1); % Index edge 23
                    % 19th B vector neighbours, G1+G2+G3
                    neighbours.Nearest(1:N-1,1:N-1,1:N-1,19) = meshInd(2:N,2:N,2:N); % Index body
                    neighbours.Nearest(1:N-1,N,1:N-1,19) = meshInd(2:N,1,2:N); % Index face 1
                    neighbours.Nearest(N,1:N-1,1:N-1,19) = meshInd(1,2:N,2:N); % Index face 2
                    neighbours.Nearest(1:N-1,1:N-1,N,19) = meshInd(2:N,2:N,1); % Index face 3
                    neighbours.Nearest(N,N,1:N-1,19) = meshInd(1,1,2:N); % Index edge 12
                    neighbours.Nearest(1:N-1,N,N,19) = meshInd(2:N,1,1); % Index edge 13
                    neighbours.Nearest(N,1:N-1,N,19) = meshInd(1,2:N,1); % Index edge 23
                    neighbours.Nearest(N,N,N,19) = meshInd(1,1,1); % Index corner
                    % 20th B vector neighbours, -(G1+G2+G3)
                    neighbours.Nearest(2:N,2:N,2:N,20) = meshInd(1:N-1,1:N-1,1:N-1); % Index body
                    neighbours.Nearest(2:N,1,2:N,20) = meshInd(1:N-1,N,1:N-1); % Index face 1
                    neighbours.Nearest(1,2:N,2:N,20) = meshInd(N,1:N-1,1:N-1); % Index face 2
                    neighbours.Nearest(2:N,2:N,1,20) = meshInd(1:N-1,1:N-1,N); % Index face 3
                    neighbours.Nearest(1,1,2:N,20) = meshInd(N,N,1:N-1); % Index edge 12
                    neighbours.Nearest(2:N,1,1,20) = meshInd(1:N-1,N,N); % Index edge 13
                    neighbours.Nearest(1,2:N,1,20) = meshInd(N,1:N-1,N); % Index edge 23
                    neighbours.Nearest(1,1,1,20) = meshInd(N,N,N); % Index corner
                    % 21st B vector neighbours, G1+G2-G3
                    neighbours.Nearest(1:N-1,1:N-1,2:N,21) = meshInd(2:N,2:N,1:N-1); % Index body
                    neighbours.Nearest(1:N-1,N,2:N,21) = meshInd(2:N,1,1:N-1); % Index face 1
                    neighbours.Nearest(N,1:N-1,2:N,21) = meshInd(1,2:N,1:N-1); % Index face 2
                    neighbours.Nearest(1:N-1,1:N-1,1,21) = meshInd(2:N,2:N,N); % Index face 3
                    neighbours.Nearest(N,N,2:N,21) = meshInd(1,1,1:N-1); % Index edge 12
                    neighbours.Nearest(1:N-1,N,1,21) = meshInd(2:N,1,N); % Index edge 13
                    neighbours.Nearest(N,1:N-1,1,21) = meshInd(1,2:N,N); % Index edge 23
                    neighbours.Nearest(N,N,1,21) = meshInd(1,1,N); % Index corner
                    % 22nd B vector neighbours, -(G1+G2-G3)
                    neighbours.Nearest(2:N,2:N,1:N-1,22) = meshInd(1:N-1,1:N-1,2:N); % Index body
                    neighbours.Nearest(2:N,1,1:N-1,22) = meshInd(1:N-1,N,2:N); % Index face 1
                    neighbours.Nearest(1,2:N,1:N-1,22) = meshInd(N,1:N-1,2:N); % Index face 2
                    neighbours.Nearest(2:N,2:N,N,22) = meshInd(1:N-1,1:N-1,1); % Index face 3
                    neighbours.Nearest(1,1,1:N-1,22) = meshInd(N,N,2:N); % Index edge 12
                    neighbours.Nearest(2:N,1,N,22) = meshInd(1:N-1,N,1); % Index edge 13
                    neighbours.Nearest(1,2:N,N,22) = meshInd(N,1:N-1,1); % Index edge 23
                    neighbours.Nearest(1,1,N,22) = meshInd(N,N,1); % Index corner
                    % 23rd B vector neighbours, G1-G2+G3
                    neighbours.Nearest(2:N,1:N-1,1:N-1,23) = meshInd(1:N-1,2:N,2:N); % Index body
                    neighbours.Nearest(2:N,N,1:N-1,23) = meshInd(1:N-1,1,2:N); % Index face 1
                    neighbours.Nearest(1,1:N-1,1:N-1,23) = meshInd(N,2:N,2:N); % Index face 2
                    neighbours.Nearest(2:N,1:N-1,N,23) = meshInd(1:N-1,2:N,1); % Index face 3
                    neighbours.Nearest(1,N,1:N-1,23) = meshInd(N,1,2:N); % Index edge 12
                    neighbours.Nearest(2:N,N,N,23) = meshInd(1:N-1,1,1); % Index edge 13
                    neighbours.Nearest(1,1:N-1,N,23) = meshInd(N,2:N,1); % Index edge 23
                    neighbours.Nearest(1,N,N,23) = meshInd(N,1,1); % Index corner
                    % 24th B vector neighbours, -(G1-G2+G3)
                    neighbours.Nearest(1:N-1,2:N,2:N,24) = meshInd(2:N,1:N-1,1:N-1); % Index body
                    neighbours.Nearest(1:N-1,1,2:N,24) = meshInd(2:N,N,1:N-1); % Index face 1
                    neighbours.Nearest(N,2:N,2:N,24) = meshInd(1,1:N-1,1:N-1); % Index face 2
                    neighbours.Nearest(1:N-1,2:N,1,24) = meshInd(2:N,1:N-1,N); % Index face 3
                    neighbours.Nearest(N,1,2:N,24) = meshInd(1,N,1:N-1); % Index edge 12
                    neighbours.Nearest(1:N-1,1,1,24) = meshInd(2:N,N,N); % Index edge 13
                    neighbours.Nearest(N,2:N,1,24) = meshInd(1,1:N-1,N); % Index edge 23
                    neighbours.Nearest(N,1,1,24) = meshInd(1,N,N); % Index corner
                    % 25th B vector neighbours, G1-G2-G3
                    neighbours.Nearest(2:N,1:N-1,2:N,25) = meshInd(1:N-1,2:N,1:N-1); % Index body
                    neighbours.Nearest(2:N,N,2:N,25) = meshInd(1:N-1,1,1:N-1); % Index face 1
                    neighbours.Nearest(1,1:N-1,2:N,25) = meshInd(N,2:N,1:N-1); % Index face 2
                    neighbours.Nearest(2:N,1:N-1,1,25) = meshInd(1:N-1,2:N,N); % Index face 3
                    neighbours.Nearest(1,N,2:N,25) = meshInd(N,1,1:N-1); % Index edge 12
                    neighbours.Nearest(2:N,N,1,25) = meshInd(1:N-1,1,N); % Index edge 13
                    neighbours.Nearest(1,1:N-1,1,25) = meshInd(N,2:N,N); % Index edge 23
                    neighbours.Nearest(1,N,1,25) = meshInd(N,1,N); % Index corner
                    % 26th B vector neighbours, -(G1-G2-G3)
                    neighbours.Nearest(1:N-1,2:N,1:N-1,26) = meshInd(2:N,1:N-1,2:N); % Index body
                    neighbours.Nearest(1:N-1,1,1:N-1,26) = meshInd(2:N,N,2:N); % Index face 1
                    neighbours.Nearest(N,2:N,1:N-1,26) = meshInd(1,1:N-1,2:N); % Index face 2
                    neighbours.Nearest(1:N-1,2:N,N,26) = meshInd(2:N,1:N-1,1); % Index face 3
                    neighbours.Nearest(N,1,1:N-1,26) = meshInd(1,N,2:N); % Index edge 12
                    neighbours.Nearest(1:N-1,1,N,26) = meshInd(2:N,N,1); % Index edge 13
                    neighbours.Nearest(N,2:N,N,26) = meshInd(1,1:N-1,1); % Index edge 23
                    neighbours.Nearest(N,1,N,26) = meshInd(1,N,1); % Index corner
                    % Remove neighbour indices with zero weight and reshape the nearest neighbour index
                    neighbours.Nearest = neighbours.Nearest(:,:,:,wbInd);
                    neighbours.Nearest = reshape(neighbours.Nearest, numQpts, numBVecs); % q,b
            end
        end
    end
end