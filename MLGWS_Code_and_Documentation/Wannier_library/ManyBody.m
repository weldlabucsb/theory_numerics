
% Calculate the hopping and many-body interaction terms for the optical lattice.
%
% PROPERTIES
% Bands: Row vector containing all the bands to calculate Wannier functions for
% Sites: Fractional coordinates of each site to calculate Wannier functions for
% W: Wannier functions, if they have been modified by a unitary transformation
% UPhase: Transformation matrix to make the Wannier functions real
% J: The hopping matrix elements given by band numbers m and n, and lattice centre j
% Ujjnn: Density-density interaction, U_00jj^mmnn
% U0jnn: Density-induced hopping within bands, U_000j^mmnn
% Ujjmn: Density-induced hopping within sites, U_00jj^mmmn
% U0jmn: Density-induced hopping: same site and band to different sites and band, U_000j^mmmn
% Uj0mn: Density-induced hopping: same band, different site to different band, same site, U_00j0^mmmn
% Uj0j0: Double hopping terms, U_j0j0^mmnn
%
% METHODS
% ManyBody(r, groups, sites): Constructor.
% CalculateWannier(manyBody, lattice, recip, bloch, wannier90, superCell, calcMode) 
% CalculateHopping(manyBody, lattice, bloch, wannier90):
% CalculateInteraction(manyBody, superCell):
%
% REFERENCES (equation numbers are referenced in the code):
% 1. D. Jaksch, C. Bruder, J.I. Cirac, C.W. Gardiner, P. Zoller, Cold bosonic atoms in optical
% lattices. Phys. Rev. Lett. 81, 3108 (1998).
%
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

classdef ManyBody
    
    properties(GetAccess = 'public', SetAccess = 'private')
        % DECLARE PROPERTIES AND DEFAULT VALUES
        Bands;
        Sites;
        SCell;
        W;
        UPhase;
        J;
        Ujjnn;
        U0jnn;
        Ujjmn;
        U0jmn;
        Uj0mn;
        Uj0j0;
    end
    
    methods
        function manyBody = ManyBody(groups, sites, superCell)
            % CLASS CONSTRUCTOR
            %
            % REQUIRED INPUT ARGUMENTS
            %    groups --> cell determining the bands belonging to each group, e.g. {1; [2 3]}
            %     sites --> manyBody.Sites, e.g. [[0; 0] [1; 0] [0; 1] [1; 1] [1; -1]];
            % superCell --> instance of the SuperCell class
            
            manyBody.Sites = sites;
            manyBody.SCell = superCell;
            % Create a row vector indexing all the bands for which to calculate Wannier functions
            bands = [];
            for n = 1 : length(groups)
                bands = [bands groups{n}]; end
            manyBody.Bands = bands;
        end
        
        function manyBody = CalculateWannier(manyBody, lattice, recip, bloch, wannier90, calcMode)
            % CALCULATE THE WANNIER FUNCTIONS IN REAL SPACE OVER A SUPERCELL
            %
            % REQUIRED INPUT ARGUMENTS
            %   lattice --> instance of the Lattice class
            %     recip --> instance of the ReciprocalLattice class
            %     bloch --> instance of the Bloch class
            % wannier90 --> instance of the Wannier90 class
            %  calcMode --> choose whether to calculate the Wannier functions using 'multiprod' or to use a 
            %               for 'loop'
            
            disp('Calculating real space Wannier data...');
            tic;
            superCell = manyBody.SCell;
            % Convert sites to cartesian coordinates
            sites = transpose(lattice.R * manyBody.Sites); 
            % Pick out the states according to manyBody.Bands
            state = bloch.State(:, :, manyBody.Bands);
            % Apply the transformation matrices wannier90.U to the Bloch states, Eq. 59 of Ref. 1
            state = multiprod(permute(state, [1 3 2]), wannier90.U);
            state = permute(state, [1 3 2]);
            % Set variables
            numBands = size(state, 3); numSites = size(sites, 1);
            numXpts = superCell.NumPts; dim = size(lattice.R, 1);
            % Calculate the Wannier functions across the superCell
            switch calcMode
                case 'multiprod'
                    % The periodic Bloch function given by q, x and n
                    u = multiprod(permute(state, [2 1 3]), exp(1i * transpose(recip.Gk) * superCell.Mesh)) ...
                        / sqrt(lattice.Vol);
                    eQX = exp(1i * transpose(bloch.Q.Mesh) * superCell.Mesh);
                    % The Bloch function given by q, x and n
                    psi = u .* eQX(:, :, ones(1, numBands)); clear eQX
                    manyBody.W = multiprod(exp(-1i * sites * bloch.Q.Mesh), psi);
                    manyBody.W = permute(manyBody.W, [2 3 1]);
                case 'loop'
                    % Using a for loop in the following manner is slower but less memory intensive than the
                    % above, and is therefore required in 3D when numKpts, numQpts, and numXpts are all
                    % large.
                    manyBody.W = zeros(numSites, numBands, numXpts);
                    for x = 1 : numXpts
                        if x / 100 == round(x / 100)
                            disp(['x = ' num2str(x) '/' num2str(numXpts)]); end
                        % The periodic Bloch function given by q, x and n
                        u = permute(multiprod(exp(1i * transpose(superCell.Mesh(:, x)) * recip.Gk), ...
                            state) / sqrt(lattice.Vol), [2 3 1]);
                        eQX = exp(1i * transpose(bloch.Q.Mesh) * superCell.Mesh(:, x));
                        % The Bloch function given by q, x and n
                        psi = u .* eQX(:, ones(1, numBands));
                        manyBody.W(:, :, x) = exp(-1i * sites * bloch.Q.Mesh) * psi;
                    end
                    manyBody.W = permute(manyBody.W, [3 2 1]);
            end
            % The Wannier function given by position x, band number n, and lattice centre j
            manyBody.W = manyBody.W * bloch.Q.dVol * lattice.Vol / ((2 * pi) ^ dim);
            % Make the Wannier functions real
            Wphase = exp(-1i * angle(squeeze(max(manyBody.W)))); % n,j
            if size(manyBody.W, 2) == 1
                Wphase = transpose(Wphase); end
            manyBody.W = manyBody.W .* permute(Wphase(:, :, ones(1, numXpts)), [3 1 2]);
            manyBody.W = real(manyBody.W);
            manyBody.UPhase = diag(Wphase(:,1));
            toc;
        end
        
        function manyBody = CalculateHopping(manyBody, lattice, bloch, wannier90)
            % CALCULATE THE HOPPING MATRIX ELEMENTS
            %
            % REQUIRED INPUT ARGUMENTS
            %   lattice --> instance of the Lattice class
            %     bloch --> instance of the Bloch class
            % wannier90 --> instance of the Wannier90 class
            
            disp('Calculating hopping matrix elements...');
            tic;
            % Convert sites to cartesian coordinates
            sites = transpose(lattice.R * manyBody.Sites); 
            % Pick out the energies according to manyBody.Bands
            energy = bloch.Energy(manyBody.Bands, :);
            % Set variables
            numQpts = size(bloch.State, 2); numBands = length(manyBody.Bands); dim = size(lattice.R, 1);
            Dq = zeros(numBands, numBands, numQpts);
            for q = 1 : numQpts
                Dq(:, :, q) = diag(energy(:, q)); end
            UDU = permute(multiprod(conj(permute(wannier90.U, [2 1 3])), multiprod(Dq, wannier90.U)), [3 1 2]); % q,m,n
            % The hopping matrix elements given by band numbers m and n, and lattice centre j
            manyBody.J = permute(multiprod(exp(1i * sites * bloch.Q.Mesh), UDU) * ...
                bloch.Q.dVol * lattice.Vol / ((2 * pi) ^ dim), [2 3 1]); % m,n,j
            toc;
        end
        
        function manyBody = CalculateInteraction(manyBody)
            % CALCULATE THE INTERACTION MATRIX ELEMENTS
            
            disp('Calculating many-body interaction matrix elements...');
            tic;
            % Set variables
            superCell = manyBody.SCell;
            numBands = length(manyBody.Bands); numSites = size(manyBody.Sites, 2);
            manyBody.Ujjnn = zeros(numBands, numBands, numSites);
            manyBody.U0jnn = zeros(numBands, numBands, numSites);
            manyBody.Ujjmn = zeros(numBands, numBands, numSites);
            manyBody.U0jmn = zeros(numBands, numBands, numSites);
            manyBody.Uj0mn = zeros(numBands, numBands, numSites);
            manyBody.Uj0j0 = zeros(numBands, numBands, numSites);
            % Calculate interaction matrix elements
            for j = 1 : numSites
                % Density-density interaction, U_00jj^mmnn
                manyBody.Ujjnn(:,:,j) = superCell.dVol * transpose(manyBody.W(:,:,1) .* manyBody.W(:,:,1)) * ...
                    (manyBody.W(:,:,j) .* manyBody.W(:,:,j));
                % Density-induced hopping within bands, U_000j^mmnn
                manyBody.U0jnn(:,:,j) = superCell.dVol * transpose(manyBody.W(:,:,1) .* manyBody.W(:,:,1)) * ...
                    (manyBody.W(:,:,1) .* manyBody.W(:,:,j));
                % Density-induced hopping within sites, U_00jj^mmmn
                manyBody.Ujjmn(:,:,j) = superCell.dVol * transpose(manyBody.W(:,:,1) .* manyBody.W(:,:,1) .* ...
                    manyBody.W(:,:,j)) * manyBody.W(:,:,j);
                % Density-induced hopping: same site and band to different sites and band, U_000j^mmmn
                manyBody.U0jmn(:,:,j) = superCell.dVol * transpose(manyBody.W(:,:,1) .* manyBody.W(:,:,1) .* ...
                    manyBody.W(:,:,1)) * manyBody.W(:,:,j);
                % Density-induced hopping: same band, different site to different band, same site, U_00j0^mmmn
                manyBody.Uj0mn(:,:,j) = superCell.dVol * transpose(manyBody.W(:,:,1) .* manyBody.W(:,:,1) .* ...
                    manyBody.W(:,:,j)) * manyBody.W(:,:,1);
                % Double hopping terms, U_j0j0^mmnn
                manyBody.Uj0j0(:,:,j) = superCell.dVol * transpose(manyBody.W(:,:,j) .* manyBody.W(:,:,1)) * ...
                    (manyBody.W(:,:,j) .* manyBody.W(:,:,1));
            end
            toc;
        end
        
        
        function Plot(manyBody, site, isoSurf, points)
            % PLOT THE WANNIER FUNCTIONS
            %
            % REQUIRED INPUT ARGUMENTS
            % site --> the site at which to plot the Wannier functions
            %
            % OPTIONAL INPUT ARGUMENTS
            % isosurf --> value of the iso-surface for a 3D plot
            
            if nargin < 2
                isoSurf = 1; end
            superCell = manyBody.SCell;
            numBands = size(manyBody.W, 2); numSites = size(manyBody.W, 3);
            switch superCell.Dimension
                case '1D'
                    for n = 1 : numBands
                        figure;
                        % hold on
                        plot(superCell.MeshX, manyBody.W(:, n, site), 'b-');
                        axis([min(superCell.MeshX) max(superCell.MeshX) min(manyBody.W(:, n, site)) ...
                            max(manyBody.W(:, n, site))]);
                        xlabel('x/\lambda'); ylabel('W*\lambda^(1/2)');
                        title(['Wannier function ' num2str(n)]);
                    end
                case '2D'
                    Wann = reshape(manyBody.W, superCell.NumXpts, superCell.NumYpts, numBands, numSites);
                    for n = 1 : numBands
                        figure;
%                         v = linspace(min(min(Wann(:, :, n, site))), max(max(Wann(:, :, n, site))), 20);
%                         contourf(superCell.MeshX, superCell.MeshY, Wann(:, :, n, site), v);
                        surf(superCell.MeshX,superCell.MeshY,Wann(:,:,n,site));
                        axis([min(superCell.MeshX) max(superCell.MeshX) min(superCell.MeshY) ...
                            max(superCell.MeshY)]);
%                         surf(superCell.MeshX, superCell.MeshY, Wann(:, :, n, site));
%                         axis([min(superCell.MeshX) max(superCell.MeshX) min(superCell.MeshY) ...
%                             max(superCell.MeshY) -Inf Inf]);
                        xlabel('x/\lambda'); ylabel('y/\lambda'); zlabel('W*\lambda');
                        title(['Wannier function ' num2str(n)]);
                        colorbar;
                        caxis([min(min(Wann(:, :, n, site))),max(max(Wann(:, :, n, site)))]);
                        if nargin > 3
                            hold on
                            plot(points(1,:), points(2,:), 'wo', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
                            hold off
                        end
                    end
                    
                case '3D'
                    Wann = reshape(manyBody.W, superCell.NumXpts, superCell.NumYpts, superCell.NumZpts, ...
                        numBands, numSites);
                    for n = 1 : numBands
                        figure;
                        fv = patch(isosurface(superCell.MeshX, superCell.MeshY, superCell.MeshZ, ...
                            abs(Wann(:, :, :, n, site)), isoSurf));
                        axis([min(superCell.MeshX) max(superCell.MeshX) min(superCell.MeshY) ...
                            max(superCell.MeshY) min(superCell.MeshZ) max(superCell.MeshZ)]);
                        xlabel('x/\lambda'); ylabel('y/\lambda'); zlabel('z/\lambda');
                        title(['Wannier function W*\lambda^(3/2)' num2str(n)]);
                        set(fv, 'FaceColor', 'red', 'EdgeColor', 'none')
                        view([-65, 20])
                        camlight left;
                        set(gcf, 'Renderer', 'zbuffer'); lighting phong
                        if nargin > 3
                            hold on
                            plot3(points(1,:), points(2,:), points(3,:), 'wo', 'MarkerFaceColor', 'k', ...
                                'MarkerSize', 5);
                            hold off
                        end
                    end
            end
        end
    end
end