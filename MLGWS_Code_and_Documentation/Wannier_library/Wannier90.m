
% Properties and methods to compute and reduce the spread of a set of Wannier functions
% 
% PROPERTIES
% Groups: The band indices grouped into the respective composite groups
% GroupsAsc: A cell equivalent to groups but with the band indices starting from 1
% Mmn: The overlap matrices for neighbouring Bloch states
% Omega: The total spread of each individual Wannier function
% OmegaI: The invariant part of the spread of each individual Wannier function
% OmegaD: The band-diagonal part of the spread of each individual Wannier function
% Centre: The position of the centre of each individual Wannier function
% U: The update matrices to transform the initial Bloch states to give maximally
%    localised Wannier functions
% GroupOmega: The total spread of each composite group of Wannier functions
% GroupOmegaI: The invariant part of the spread of each composite group of Wannier
%              functions
% GroupOmegaD: The band-diagonal part of the spread of each composite group of Wannier
%              functions
% GroupOmegaOD: The band-off-diagonal part of the spread of each composite group of
%              Wannier functions
% 
% METHODS
% Wannier90(state, meshInd, recip, neighbours): Constructor.
% Isolated(wannier90, dimension, mesh, neighbours, epsilon, maxIter)
% Composite(wannier90, mesh, neighbours, epsilon, maxIter)
% Disentangle(wannier90, neighbours, alpha, maxIter, randomise)
% SingleBandProperties(wannier90, neighbours)
% GroupProperties(wannier90, neighbours)
% SetGroups(wannier90, groups)
% Transform(wannier90, neighbours, UTrans)
% 
% REFERENCES (equation numbers are referenced in the code):
% 1. N. Marzari and D. Vanderbilt. Maximally localized generalized Wannier functions for composite energy
% bands. Phys. Rev. B 56(20), 12847-12865 Nov 1997.
% 2. I. Souza, N. Marzari and D. Vanderbilt. Maximally localized Wannier functions for entangled energy bands.
% Phys. Rev. B 65(3), 035109 Dec 2001.
% 3. A. Mostofi, J. Yates, Y-S Lee, I. Souza, D. Vanderbilt and N. Marzari.  wannier90: A tool for obtaining
% maximally-localised Wannier functions. Comp. Phys. Comm 178(9) 685-699 May 2008
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch. 
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

classdef Wannier90
    
    properties(GetAccess = 'public', SetAccess = 'public')
        % DECLARE PROPERTIES AND DEFAULT VALUES
        Groups;
        GroupsAsc;
        Mmn;
        Omega;
        OmegaI;
        OmegaD;
        Centre;
        U;
        GroupOmega;
        GroupOmegaI;
        GroupOmegaD;
        GroupOmegaOD;
    end
    
    methods
        
        function wannier90 = Wannier90(state, meshInd, recip, neighbours)
            % CLASS CONSTRUCTOR
            %
            % REQUIRED INPUT ARGUMENTS
            %      state --> Fourier coefficients of the Bloch states, labelled by k
            %    meshInd --> indices of the Brillouin zone mesh points
            %      recip --> instance of the ReciprocalLattice class
            % neighbours --> instance of the Neighbours class
            
            dim = size(neighbours.B, 1); numQpts = size(state, 2); numBands = size(state, 3);
            wannier90.Groups = cell(numBands, 1);
            for n = 1 : numBands
                wannier90.Groups{n} = n; end
            wannier90.GroupsAsc = wannier90.Groups;
            [wannier90.Mmn] = OverlapsComposite(state, meshInd, recip, neighbours);
            wannier90.Omega = zeros(1, numBands);
            wannier90.OmegaI = zeros(1, numBands);
            wannier90.OmegaD = zeros(1, numBands);
            wannier90.Centre = zeros(dim, numBands);
            wannier90.U = eye(numBands);
            wannier90.U = wannier90.U(:, :, ones(1, numQpts));
        end
        
        function wannier90 = Isolated(wannier90, dimension, mesh, neighbours, epsilon, maxIter)
            % IMPLEMENTS THE ISOLATED BAND WANNIER90 ALGORITHM ON THE INSTANCE
            %
            % REQUIRED INPUT ARGUMENTS
            %  dimension --> dimension of the system, e.g. '1D'
            %       mesh --> a mesh of Brillouin zone points
            % neighbours --> instance of the Neighbours class
            %    epsilon --> a parameter used in the Wannier90 algorithm. Set between 0 and 1
            %    maxIter --> the maximum number of Wannier90 iterations
            
            disp('Minimising the Wannier spread of the Bloch bands, treating each as isolated...');
            tic;
            numBvecs = size(neighbours.B, 2);
            numBands = size(wannier90.Mmn, 1); numQpts = size(wannier90.Mmn, 3);
            UTemp = eye(numBands); UTemp = UTemp(:, :, ones(1, numQpts));
            % Use the steepest descent algorithm to minimise omegaD
            for n = 1 : numBands
                disp(['Minimising omegaD for Bloch band ' num2str(n)]);
                % If the Wannier states have inversion symmetry, omegaD = 0, and the phases can be set manually
                % (see Sec. IV.C.3 of Ref 1)
                [UInv, MmnTemp] = Wannier90Inversion(dimension, neighbours, squeeze(wannier90.Mmn(n,n,:,:)));
                % Run the Wannier90 algorithm outlined in Ref. 1 for an isolated band
                [UIso, omega, omegaI, omegaD, centre] = Wannier90Isolated(MmnTemp, neighbours, epsilon, ...
                    maxIter, mesh);
                % Save the variables for each band
                UTemp(n,n,:) = UIso .* UInv;
                wannier90.Omega(n) = omega(end);
                wannier90.OmegaI(n) = omegaI(end);
                wannier90.OmegaD(n) = omegaD(end);
                wannier90.Centre(:,n) = centre;
            end
            wannier90.U = multiprod(wannier90.U, UTemp);
            % Update the overlaps between nn. mesh points, Eq. 61
            wannier90.Mmn = multiprod(conj(permute(UTemp(:, :, :, ones(1, numBvecs)), [2 1 3 4])), ...
                wannier90.Mmn);
            wannier90.Mmn = multiprod(wannier90.Mmn, reshape(UTemp(:,:,neighbours.Nearest), numBands, ...
                numBands, numQpts, numBvecs));
            toc;
        end
        
        function wannier90 = ParallelTransport(wannier90, mesh, neighbours)
            % IMPLEMENTS THE PARALLEL TRANSPORT ALGORITHM ON THE INSTANCE (1D ONLY)
            %
            % REQUIRED INPUT ARGUMENTS
            %       mesh --> a mesh of Brillouin zone points
            % neighbours --> instance of the Neighbours class
            
            if size(mesh, 1) ~= 1
                disp('ParallelTransport may only be used in one dimension');
                return;
            end
            disp('Using parallel transport to minimise the Wannier spread of each composite group...');
            tic;
            groups = wannier90.GroupsAsc; dim = size(neighbours.B, 1);
            numBands = size(wannier90.Mmn, 1); numQpts = size(wannier90.Mmn, 3);
            UComp = eye(numBands); UComp = UComp(:, :, ones(1, numQpts));
            % Use the steepest descent algorithm to minimise omegaOD
            for n = 1 : length(groups)
                if length(groups{n}) > 1
                    disp(['Band group ' num2str(n) ': minimising the gauge-dependent part of the spread \tilde{omega} via steepest-descent...']);
                    % Run the Wannier90 algorithm for the composite group
                    [UTemp, omega, omegaI, omegaD, omegaOD, centre, MmnTemp] = Wannier90ParallelTransport(...
                        wannier90.Mmn(groups{n}, groups{n}, :, :), neighbours, mesh);
                    % Save the variables for each group
                    wannier90.Mmn(groups{n}, groups{n}, :, :) = MmnTemp;
                    UComp(groups{n}, groups{n}, :) = UTemp;
                    wannier90.GroupOmega(n) = omega(end);
                    wannier90.GroupOmegaI(n) = omegaI(end);
                    wannier90.GroupOmegaD(n) = omegaD(end);
                    wannier90.GroupOmegaOD(n) = omegaOD(end);
                else % Isolated bands should already be minimised using the isolated version above
                    disp(['Band group ' num2str(n) ': isolated band']);
                    % The phases between nearest-neighbour q-points
                    phase = imag(log(squeeze(wannier90.Mmn(n,n,:,:))));
                    % Calculate the position of the Wannier centre, Eq. 31 of Ref. 1
                    centre = -sum(repmat(neighbours.Weight .* sum(phase, 1), dim, 1) .* neighbours.B, 2) / ...
                        numQpts;
                    % Calculate the diagonal part of the Wannier spread, Eq. 36 of Ref. 1
                    Brdash = transpose(centre) * neighbours.B;
                    qn = phase + Brdash(ones(1, numQpts), :);
                    wannier90.GroupOmegaD(n) = sum(((-qn) .^ 2) * transpose(neighbours.Weight)) / numQpts;
                    % Calculate the invariant part of the Wannier spread, Eq. 34 of Ref. 1
                    wannier90.GroupOmegaI(n) = neighbours.Weight * transpose(1 - sum(abs(squeeze( ...
                        wannier90.Mmn(n,n,:,:))) .^ 2, 1) / numQpts);
                    wannier90.GroupOmega(n) = wannier90.OmegaI(n) + wannier90.OmegaD(n);
                end
            end
            wannier90.U = multiprod(wannier90.U, UComp);
            wannier90 = wannier90.SingleBandProperties(neighbours);
            toc;
        end
        
        function wannier90 = Composite(wannier90, mesh, neighbours, epsilon, maxIter)
            % IMPLEMENTS THE COMPOSITE BAND WANNIER90 ALGORITHM ON THE INSTANCE
            %
            % REQUIRED INPUT ARGUMENTS
            %       mesh --> a mesh of Brillouin zone points
            % neighbours --> instance of the Neighbours class
            %    epsilon --> a parameter used in the Wannier90 algorithm. Set between 0 and 1
            %    maxIter --> the maximum number of Wannier90 iterations
            
            disp('Minimising the Wannier spread of each composite group...');
            tic;
            groups = wannier90.GroupsAsc; dim = size(neighbours.B, 1);
            numBands = size(wannier90.Mmn, 1); numQpts = size(wannier90.Mmn, 3);
            UComp = eye(numBands); UComp = UComp(:, :, ones(1, numQpts));
            % Use the steepest descent algorithm to minimise omegaOD
            for n = 1 : length(groups)
                if length(groups{n}) > 1
                    disp(['Band group ' num2str(n) ': minimising the gauge-dependent part of the spread \tilde{omega} via steepest-descent...']);
                    % Run the Wannier90 algorithm for the composite group
                    [UTemp, omega, omegaI, omegaD, omegaOD, centre, MmnTemp] = Wannier90Composite(...
                        wannier90.Mmn(groups{n}, groups{n}, :, :), neighbours, epsilon, maxIter, mesh);
                    % Save the variables for each group
                    wannier90.Mmn(groups{n}, groups{n}, :, :) = MmnTemp;
                    UComp(groups{n}, groups{n}, :) = UTemp;
                    wannier90.GroupOmega(n) = omega(end);
                    wannier90.GroupOmegaI(n) = omegaI(end);
                    wannier90.GroupOmegaD(n) = omegaD(end);
                    wannier90.GroupOmegaOD(n) = omegaOD(end);
                else % Isolated bands should already be minimised using the isolated version above
                    disp(['Band group ' num2str(n) ': isolated band']);
                    % The phases between nearest-neighbour q-points
                    phase = imag(log(squeeze(wannier90.Mmn(n,n,:,:))));
                    % Calculate the position of the Wannier centre, Eq. 31 of Ref. 1
                    centre = -sum(repmat(neighbours.Weight .* sum(phase, 1), dim, 1) .* neighbours.B, 2) / ...
                        numQpts;
                    % Calculate the diagonal part of the Wannier spread, Eq. 36 of Ref. 1
                    Brdash = transpose(centre) * neighbours.B;
                    qn = phase + Brdash(ones(1, numQpts), :);
                    wannier90.GroupOmegaD(n) = sum(((-qn) .^ 2) * transpose(neighbours.Weight)) / numQpts;
                    % Calculate the invariant part of the Wannier spread, Eq. 34 of Ref. 1
                    wannier90.GroupOmegaI(n) = neighbours.Weight * transpose(1 - sum(abs(squeeze( ...
                        wannier90.Mmn(n,n,:,:))) .^ 2, 1) / numQpts);
                    wannier90.GroupOmega(n) = wannier90.OmegaI(n) + wannier90.OmegaD(n);
                end
            end
            wannier90.U = multiprod(wannier90.U, UComp);
            wannier90 = wannier90.SingleBandProperties(neighbours);
            toc;
        end
        
        function wannier90 = Disentangle(wannier90, neighbours, alpha, maxIter, randomise)
            % IMPLEMENTS THE DISENTANGLING ALGORITHM ON THE INSTANCE
            %
            % REQUIRED INPUT ARGUMENTS
            % neighbours --> instance of the Neighbours class
            %      alpha --> a parameter used in the algorithm. Set between 0 and 1
            %    maxIter --> the maximum number of iterations
            %  randomise --> set 'True' to randomise the initial configuration
            
            disp('Using the Souza et al. algorithm to reduce the contribution to omega_OD from each group of bands...');
            tic;
            groups = wannier90.GroupsAsc; dim = size(neighbours.B, 1);
            numBands = size(wannier90.Mmn, 1); numQpts = size(wannier90.Mmn, 3);
            UComp = eye(numBands); UComp = UComp(:, :, ones(1, numQpts));
            % Use the disentangling algorithm to minimise omegaI for each band in the group
            for n = 1 : length(groups)
                if length(groups{n}) > 1 % Disentangle each band in the group in turn
                    disp(['Band group ' num2str(n) ': optimally extacting a single band...']);
                    % Run the wannier90 disentangling algorithm to minimise omegaI for each band
                    [UTemp, omega, omegaI, omegaD, omegaOD, centre, MmnTemp] = Wannier90Disentangle( ...
                        wannier90.Mmn(groups{n}, groups{n}, :, :), neighbours, alpha, maxIter, randomise);
                    % Save the variables for each group
                    wannier90.Mmn(groups{n}, groups{n}, :, :) = MmnTemp;
                    UComp(groups{n}, groups{n}, :) = UTemp;
                    wannier90.GroupOmega(n) = omega(end);
                    wannier90.GroupOmegaI(n) = omegaI(end);
                    wannier90.GroupOmegaD(n) = omegaD(end);
                    wannier90.GroupOmegaOD(n) = omegaOD(end);
                else % Isolated bands don't need disentangling
                    disp(['Band group ' num2str(n) ': isolated band']);
                    % The phases between nearest-neighbour q-points
                    phase = imag(log(squeeze(wannier90.Mmn(n,n,:,:))));
                    % Calculate the position of the Wannier centre, Eq. 31 of Ref. 1
                    centre = -sum(repmat(neighbours.Weight .* sum(phase, 1), dim, 1) .* neighbours.B, 2) / ...
                        numQpts;
                    % Calculate the diagonal part of the Wannier spread, Eq. 36 of Ref. 1
                    Brdash = transpose(centre) * neighbours.B;
                    qn = phase + Brdash(ones(1, numQpts), :);
                    wannier90.GroupOmegaD(n) = sum(((-qn) .^ 2) * transpose(neighbours.Weight)) / numQpts;
                    % Calculate the invariant part of the Wannier spread, Eq. 34 of Ref. 1
                    wannier90.GroupOmegaI(n) = neighbours.Weight * transpose(1 - sum(abs(squeeze( ...
                        wannier90.Mmn(n,n,:,:))) .^ 2, 1) / numQpts);
                    wannier90.GroupOmega(n) = wannier90.OmegaI(n) + wannier90.OmegaD(n);
                end
            end
            wannier90.U = multiprod(wannier90.U, UComp);
            wannier90 = wannier90.SingleBandProperties(neighbours);
            toc;
        end
        
        function wannier90 = SingleBandProperties(wannier90, neighbours)
            % SET THE SINGLE BAND PROPERTIES OF THE INSTANCE
            %
            % REQUIRED INPUT ARGUMENTS
            % neighbours --> instance of the Neighbours class
            
            dim = size(neighbours.B, 1);
            numBands = size(wannier90.Mmn, 1); numQpts = size(wannier90.Mmn, 3);
            for n = 1 : numBands
                % The phases between nearest-neighbour q-points
                phase = imag(log(squeeze(wannier90.Mmn(n,n,:,:))));
                % Calculate the position of the Wannier centre, Eq. 31 of Ref. 1
                wannier90.Centre(:,n) = -sum(repmat(neighbours.Weight .* sum(phase, 1), dim, 1) .* ...
                    neighbours.B, 2) / numQpts;
                % Calculate the diagonal part of the Wannier spread, Eq. 36 of Ref. 1
                Brdash = transpose(wannier90.Centre(:,n)) * neighbours.B;
                qn = phase + Brdash(ones(1, numQpts), :);
                wannier90.OmegaD(n) = sum(((-qn) .^ 2) * transpose(neighbours.Weight)) / numQpts;
                % Calculate the invariant part of the Wannier spread, Eq. 34 of Ref. 1
                wannier90.OmegaI(n) = neighbours.Weight * transpose(1 - sum(abs( ...
                    squeeze(wannier90.Mmn(n,n,:,:))) .^ 2, 1) / numQpts);
            end
            % The total Wannier spread for each band, Eq. 13 of Ref. 1
            wannier90.Omega = wannier90.OmegaD + wannier90.OmegaI;
        end
        
        function wannier90 = GroupProperties(wannier90, neighbours)
            % SET THE GROUP PROPERTIES OF THE INSTANCE
            %
            % REQUIRED INPUT ARGUMENTS
            % neighbours --> instance of the Neighbours class
            
            dim = size(neighbours.B, 1);
            numQpts = size(wannier90.Mmn, 3);
            numBvecs = size(wannier90.Mmn, 4);
            groups = wannier90.GroupsAsc;
            
            for n = 1 : length(groups)
                numBands = length(groups{n});
                if numBands > 1
                    MmnTemp = wannier90.Mmn(groups{n},groups{n},:,:);
                    % Update the overlaps between nn. mesh points (within bands)
                    Mnn = zeros(numBands, numBands, numQpts, numBvecs);
                    for m = 1 : numBands
                        Mnn(1,m,:,:) = MmnTemp(m,m,:,:); end
                    Mnn = Mnn(ones(1, numBands),:,:,:); % n,n,q,b
                    % Update the phases between nn. mesh points, q,b,n,n
                    phase = imag(log(permute(Mnn, [3 4 2 1])));
                    phase = permute(phase(:,:,:,1), [3 2 1]); % n,b,q
                    % Weight the phases in each direction, n,b
                    wbphase = neighbours.Weight(ones(1, numBands), :) .* sum(phase, 3);
                    % Recalculate the Wannier centre, Eq. 31
                    centre = permute(wbphase(:, :, ones(1, dim)), [3,2,1]) .* neighbours.B(:, :, ...
                        ones(1, numBands)); % dim,b,n
                    centre = -sum(permute(centre, [1 3 2]), 3) / numQpts; % dim,n
                    Brdash = transpose(centre) * neighbours.B; % n,b
                    qn = phase + Brdash(:, :, ones(1, numQpts)); % Recalculate qn, Eq. 47, n,b,q
                    qn = permute(qn, [1 3 2]); % n,q,b
                    % Recalculate omegaD, Eq. 36
                    wannier90.GroupOmegaD(n) = neighbours.Weight * sum(sum(permute(qn .^ 2, [3 1 2]), 3), ...
                        2) / numQpts;
                    % Absolute value squared of the overlaps
                    Mmn2 = abs(permute(MmnTemp, [4 1 2 3])) .^ 2;
                    % Absolute value squared of the overlaps (within bands)
                    Mnn2 = abs(permute(Mnn(1,:,:,:), [4 2 3 1])) .^ 2; 
                    % Recalculate omegaOD, Eq. 35
                    wannier90.GroupOmegaOD(n) = neighbours.Weight * (sum(sum(sum(Mmn2, 4), 3), 2) - ...
                        sum(sum(Mnn2, 3), 2)) / numQpts;
                    % Calculate the invariant part of the total spread for the composite group
                    MmnTemp = permute(MmnTemp, [4,1,2,3]);
                    wannier90.GroupOmegaI(n) = neighbours.Weight * (numBands - ...
                        sum(sum(sum(conj(MmnTemp) .* MmnTemp, 4), 3), 2) / numQpts); % Eq. 34
                    % wannier90.Mmn = permute(wannier90.Mmn, [2,3,4,1]);
                else
                    % The phases between nearest-neighbour q-points
                    phase = imag(log(squeeze(wannier90.Mmn(groups{n},groups{n},:,:))));
                    % Calculate the position of the Wannier centre, Eq. 31 of Ref. 1
                    centre = -sum(repmat(neighbours.Weight .* sum(phase, 1), dim, 1) .* neighbours.B, 2) / ...
                        numQpts;
                    % Calculate the diagonal part of the Wannier spread, Eq. 36 of Ref. 1
                    Brdash = transpose(centre) * neighbours.B;
                    qn = phase + Brdash(ones(1, numQpts), :);
                    wannier90.GroupOmegaD(n) = sum(((-qn) .^ 2) * transpose(neighbours.Weight)) / numQpts;
                    % Calculate the invariant part of the Wannier spread, Eq. 34 of Ref. 1
                    wannier90.GroupOmegaI(n) = neighbours.Weight * transpose(1 - sum(abs(squeeze( ...
                        wannier90.Mmn(groups{n},groups{n},:,:))) .^ 2, 1) / numQpts);
                end
            end
            wannier90.GroupOmega = wannier90.GroupOmegaI + wannier90.GroupOmegaD + wannier90.GroupOmegaOD;
        end
        
        function wannier90 = SetGroups(wannier90, groups)
            % SET THE GROUPS OF THE INSTANCE
            %
            % REQUIRED INPUT ARGUMENTS
            % groups --> groups together bands that are degenerate, e.g. {1,[2 3]}
            
            numGroups = length(groups);
            % Create a row vector, wband_ind, indexing all the bands for which to calculate Wannier functions
            groupBands = [];
            for n = 1 : length(groups)
                groupBands = [groupBands groups{n}]; end
            % Save the Groups property
            wannier90.Groups = groups;
            % Create a cell equivalent to groups but with the band indices starting from 1.
            wannier90.GroupsAsc = cell(length(groups), 1); bandsInd = 1 : length(groupBands);
            for n = 1 : length(groups)
                wannier90.GroupsAsc{n} = bandsInd(1 : length(groups{n})); bandsInd(1 : length(groups{n})) = []; end
            % Remove excluded bands from wannier90 properties
            wannier90.Mmn = wannier90.Mmn(groupBands, groupBands, :, :);
            wannier90.Omega = wannier90.Omega(groupBands);
            wannier90.OmegaI = wannier90.OmegaI(groupBands);
            wannier90.OmegaD = wannier90.OmegaD(groupBands);
            wannier90.Centre = wannier90.Centre(:, groupBands);
            wannier90.U = wannier90.U(groupBands, groupBands, :);
            % Initialise the Group Spread properties
            wannier90.GroupOmega = zeros(1, numGroups);
            wannier90.GroupOmegaI = zeros(1, numGroups);
            wannier90.GroupOmegaD = zeros(1, numGroups);
            wannier90.GroupOmegaOD = zeros(1, numGroups);
        end
        
        function wannier90 = Transform(wannier90, neighbours, UTrans)
            % UPDATE THE OVERLAP MATRICES WITH UTRANS ANR RESET THE INSTANCE PROPERTIES
            %
            % REQUIRED INPUT ARGUMENTS
            % neighbours --> instance of the Neighbours class
            %     UTrans --> the update matrix
            
            numBands = size(wannier90.Mmn, 1); numQpts = size(wannier90.Mmn, 3);
            numBvecs = size(wannier90.Mmn, 4);
            UTrans = UTrans(:, :, ones(1, numQpts));
            % Update the transformation matrices
            wannier90.U = multiprod(wannier90.U, UTrans);
            % Update the overlaps between nn. mesh points, Eq. 61
            wannier90.Mmn = multiprod(conj(permute(UTrans(:, :, :, ones(1, numBvecs)), [2 1 3 4])), ...
                wannier90.Mmn);
            wannier90.Mmn = multiprod(wannier90.Mmn, reshape(UTrans(:,:,neighbours.Nearest), numBands, ...
                numBands, numQpts, numBvecs));
            % Update the single band properties
            wannier90 = wannier90.SingleBandProperties(neighbours);
            % Update the group properties
            wannier90 = wannier90.GroupProperties(neighbours);
        end
        
        function wannier90 = Shift(wannier90, mesh, neighbours, shiftVecs)
            % UPDATE THE OVERLAP MATRICES USING SHIFTVECS TO SHIFT THE
            % WANNIER FUNCTIONS TO DIFFERENT LATTICE CELLS
            %
            % REQUIRED INPUT ARGUMENTS
            %       mesh --> a mesh of Brillouin zone points
            % neighbours --> instance of the Neighbours class
            %  shiftVecs --> a matrix of vectors for shifting each Wannier
            %                function in the group, e.g. [-lattice.R(:,2), [0;0]];
            
            numBands = size(wannier90.Mmn, 1); numQpts = size(wannier90.Mmn, 3);
            numBvecs = size(wannier90.Mmn, 4);
            UTemp = exp(-1i * transpose(shiftVecs) * mesh);
            UTemp = UTemp(:,:,ones(1, numBands));
            UTemp = permute(UTemp, [1 3 2]);
            UTrans = eye(numBands);
            UTrans = UTrans(:, :, ones(1, numQpts));
            UTrans = UTrans .* UTemp;
            % Update the transformation matrices
            wannier90.U = multiprod(wannier90.U, UTrans);
            % Update the overlaps between nn. mesh points, Eq. 61
            wannier90.Mmn = multiprod(conj(permute(UTrans(:, :, :, ones(1, numBvecs)), [2 1 3 4])), ...
                wannier90.Mmn);
            wannier90.Mmn = multiprod(wannier90.Mmn, reshape(UTrans(:,:,neighbours.Nearest), numBands, ...
                numBands, numQpts, numBvecs));
            % Update the single band properties
            wannier90 = wannier90.SingleBandProperties(neighbours);
            % Update the group properties
            wannier90 = wannier90.GroupProperties(neighbours);
        end
    end
end