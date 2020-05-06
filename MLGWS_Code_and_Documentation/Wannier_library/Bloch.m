
% Calculates the Bloch states and energies given a lattice, potential, and reciprocal lattice mesh
%
% PROPERTIES
% State: The Fourier coefficients of the Bloch states, labelled by k
% Energy: The band energies, labelled by n,q
% Hv: The potential energy term of the Hamiltonian
% Q: an instance of the SuperCell class, giving a uniform mesh of independent q-points on which to
% calculate the Bloch states and energies.
%
% METHODS
% Bloch(lattice, recip, potential): Constructor
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

classdef Bloch
    
    % DECLARE PROPERTIES AND DEFAULT VALUES
    properties(GetAccess = 'public', SetAccess = 'public')
        State;
        Energy;
        Hv;
        Q;
    end
    
    methods
        function bloch = Bloch(lattice, recip, potential, timeReversal, N, numBands)
            % CLASS CONSTRUCTOR
            %
            % REQUIRED INPUT ARGUMENTS
            %      lattice --> an instance of the Lattice class
            %        recip --> an instance of the ReciprocalLattice class
            %    potential --> an instance of the Potential class
            % timeReversal --> choose whether to use time-reversal symmetry
            %                  'origin' : central q-point is at the origin
            %                  'shifted' : central q-point is shifted to [1/N 1/N 1/N]
            %                  'false' : do not use time-reversal symmetry
            %            N --> number of q-points in each directions
            %     numBands --> number of bands to calculate Bloch data for
            
% ---------------------------------------------------------
% Calculate the potential energy term Hv of the Hamiltonian
% ---------------------------------------------------------
            vk = potential.Vk; numKpts = recip.NumKpts; gkm = recip.Gkm;
            vk(numKpts + 1) = 0; gkm(gkm == 0) = numKpts + 1;
            bloch.Hv = vk(gkm) / sqrt(lattice.Vol);
            keCoeff = 1 / (4 * pi ^ 2); % Units are lambda = 1, E_R = 1, hbar = 1
            
% ------------------------------------------------------------------------------------------------------------
% If the central q-point is at the origin, use time-reversal symmetry to obtain Bloch data for one half of the
% mesh from the other half
% ------------------------------------------------------------------------------------------------------------
            disp('Calculating the bandstructure across a uniform mesh...');
            tic;
            if strcmp(timeReversal, 'origin')
                disp('Using time-reversal symmetry; mesh centred on origin');
                % Create a mesh of quasi-momenta points
                qTemp = SuperCell(lattice.Dimension, lattice.G, [N N N], [0 0; 0 0; 0 0], 'true', [0 0 0]);
                % Set variables
                numQpts = qTemp.NumPts;
                qk = permute(qTemp.Mesh(:, 1 : (numQpts + 1) / 2, ones(1, numKpts)), [1 3 2]);
                energy = zeros(numBands, numQpts);
                state = zeros(numKpts, numBands, numQpts);
                % Loop over q-points
                for q = 1 : (numQpts + 1) / 2
                    if (q - 1) / 10 == round((q - 1) / 10)
                        disp(['q = ' num2str(q) '/' num2str((numQpts + 1) / 2)]); end
                    % The kinetic energy term of the Hamiltonian
                    h = keCoeff * diag(sum((recip.Gk + qk(:,:,q)) .^ 2, 1));
                    % The full Hamiltonian at q-point q
                    h = h + bloch.Hv;
                    % Diagonalise the Hamiltonian
                    [vec, en] = eig(h); 
                    % Sort the energies into ascending order.
                    [en, ind] = sort(real(diag(en)), 'ascend');
                    % The band energies
                    energy(:,q) = en(1 : numBands);
                    % The Fourier coefficients state_q,k of the periodic function u_q(x)
                    state(:,:,q) = vec(:, ind(1 : numBands));
                end
                % Infer the energy and states for the other half of the mesh
                energy(:, 1 + (numQpts + 1) / 2 : numQpts) = ... 
                    energy(:, -1 + (numQpts + 1) / 2 : -1 : 1);
                state(:, :, 1 + (numQpts + 1) / 2 : numQpts) = ...
                    conj(state(numKpts : -1 : 1, :, -1 + (numQpts + 1) / 2 : -1 : 1));
                % Create a new mesh that doesn't overlap at the edges
                bloch.Q = SuperCell(lattice.Dimension, lattice.G, [N N N], [0 0; 0 0; 0 0], 'false', [0 0 0]);
                numQpts = bloch.Q.NumPts;
                % Remove energies and states for overlapping q-points
                switch lattice.Dimension
                    case '1D'
                        energy = energy(:, 1 : end - 1);
                        state = state(:, :, 1 : end - 1);
                    case '2D'
                        energy = reshape(energy, numBands, N + 1, N + 1);
                        energy = energy(:, 1 : N, 1 : N);
                        energy = reshape(energy, numBands, numQpts);
                        state = reshape(state, numKpts, numBands, N + 1, N + 1);
                        state = state(:, :, 1 : N, 1 : N);
                        state = reshape(state, numKpts, numBands, numQpts);
                    case '3D'
                        energy = reshape(energy, numBands, N + 1, N + 1, N + 1);
                        energy = energy(:, 1 : N, 1 : N, 1 : N);
                        energy = reshape(energy, numBands, numQpts);
                        state = reshape(state, numKpts, numBands, N + 1, N + 1, N + 1);
                        state = state(:, :, 1 : N, 1 : N, 1 : N);
                        state = reshape(state, numKpts, numBands, numQpts);
                end
                
% --------------------------------------------------------------------------------------------------------------
% If the central q-point is shifted from the origin by [1/N 1/N 1/N], use time-reversal symmetry to obtain Bloch
% data for one half of the mesh from the other half
% --------------------------------------------------------------------------------------------------------------
            elseif strcmp(timeReversal, 'shifted')
                disp('Using time-reversal symmetry; mesh shifted from origin');
                % Create a mesh of quasi-momenta points
                bloch.Q = SuperCell(lattice.Dimension, lattice.G, [N N N], [0 0; 0 0; 0 0], 'false', [1/N 1/N 1/N]);
                numQpts = bloch.Q.NumPts;
                % Set variables
                qk = permute(bloch.Q.Mesh(:, :, ones(1, numKpts)), [1 3 2]);
                energy = zeros(numBands, numQpts);
                state = zeros(numKpts, numBands, numQpts);
                % Loop over q-points
                for q = 1 : numQpts / 2
                    if (q - 1) / 10 == round((q - 1) / 10)
                        disp(['q = ' num2str(q) '/' num2str(numQpts / 2)]); end
                    % The kinetic energy term of the Hamiltonian
                    h = keCoeff * diag(sum((recip.Gk + qk(:,:,q)) .^ 2, 1));
                    % The full Hamiltonian at q-point q
                    h = h + bloch.Hv;
                    % Diagonalise the Hamiltonian
                    [vec, en] = eig(h); 
                    % Sort the energies into ascending order.
                    [en, ind] = sort(real(diag(en)), 'ascend');
                    % The band energies
                    energy(:,q) = en(1 : numBands);
                    % The Fourier coefficients state_q,k of the periodic function u_q(x)
                    state(:,:,q) = vec(:, ind(1 : numBands));
                end
                % Infer the energy and states for the other half of the mesh
                energy(:, 1 + (numQpts) / 2 : numQpts) = energy(:, numQpts / 2: -1 : 1);
                state(:, :, 1 + (numQpts) / 2 : numQpts) = ... 
                    conj(state(numKpts : -1 : 1, :, numQpts / 2: -1 : 1));
                
% ------------------------------------------------------------------------------------------
% Otherwise, calculate the Bloch data at all mesh-points (cannot use time-reversal symmetry)
% ------------------------------------------------------------------------------------------
            else
                disp('Not using time-reversal symmetry');
                % Create a mesh of quasi-momenta points
                bloch.Q = SuperCell(lattice.Dimension, lattice.G, [N N N], [0 0; 0 0; 0 0], 'false', [1/N 1/N 1/N]);
                numQpts = bloch.Q.NumPts;
                % Set variables
                qk = permute(bloch.Q.Mesh(:, :, ones(1, numKpts)), [1 3 2]);
                energy = zeros(numBands, numQpts);
                state = zeros(numKpts, numBands, numQpts);
                % Loop over q-points
                for q = 1 : num_qpts
                    if (q - 1) / 10 == round((q - 1) / 10)
                        disp(['q = ' num2str(q) '/' num2str(numQpts)]); end
                    % The kinetic energy term of the Hamiltonian
                    h = keCoeff * diag(sum((recip.Gk + qk(:,:,q)) .^ 2, 1));
                    % The full Hamiltonian at q-point q
                    h = h + bloch.Hv;
                    % Diagonalise the Hamiltonian
                    [vec, en] = eig(h); 
                    % Sort the energies into ascending order.
                    [en, ind] = sort(real(diag(en)), 'ascend');
                    % The band energies
                    energy(:,q) = en(1 : numBands);
                    % The Fourier coefficients state_q,k of the periodic function u_q(x)
                    state(:,:,q) = vec(:, ind(1 : numBands));
                end
            end
            % Set the bloch properties for the state and energy
            state = permute(state, [1 3 2]);
            state(abs(state) < 1e-8) = state(abs(state) < 1e-8) + 1e-30i;
            bloch.State = state;
            bloch.Energy = energy;
            toc;
        end
    end
end