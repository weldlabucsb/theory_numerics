
% Defines an optical lattice potential by its Fourier coefficients and
% corresponding HKL matrix.
%
% PROPERTIES
% Type: Type of lattice potential
% Vi: A vector holding the values of the potential coefficients, ordered by hkl
% HKL: The reciprocal lattice points corresponding to potential coefficients, given in units of
% the reciprocal lattice vectors
% V0: The depth of the lattice potential
% Vk: Fourier coefficients taken over a reciprocal lattice mesh {Gk}
%
% METHODS
% Potential(vi, hkl, v0): Constructor
% Coefficients(potential, lattice, G, recip): Calculates the Fourier coefficients of the potential
% given a lattice object, reciprocal lattice vectors and a reciprocal lattice object.
% Plot(potential, lattice, recip, superCell, isoSurf): Plot the potential
%
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

classdef Potential
    
    properties(GetAccess = 'public', SetAccess = 'private')
        % DECLARE PROPERTIES AND DEFAULT VALUES
        Vi;
        HKL;
        V0;
        Vk;
    end
    
    methods
        function potential = Potential(vi, hkl, v0)
            % CLASS CONSTRUCTOR
            %
            % REQUIRED INPUT ARGUMENTS
            %  vi --> potential.Vi
            %  v0 --> potential.V0
            % hkl --> potential.HKL
            
            potential.Vi = vi;
            potential.HKL = hkl;
            potential.V0 = v0;
        end
        
        function potential = Coefficients(potential, lattice, G, recip)
            % SET THE COEFFICIENTS
            %
            % REQUIRED INPUT ARGUMENTS
            %   lattice --> instance of the Lattice class
            %     recip --> instance of the ReciprocalLattice class
            
            kCentre = (recip.NumKpts + 1) / 2;
            % A vector holding the values of the potential coefficients, in
            % order of reciprocal lattice mesh index
            vk = zeros(1, recip.NumKpts);
            % Coordinates of the potential coefficients
            gPot = G * potential.HKL;
            
            for n = 1 : size(gPot, 2)
                % Find the mesh index of the coefficient coordinate
                vkInd = sum((recip.Gk - gPot(:, n*ones(1, recip.NumKpts))) .^ 2, 1);
                % Set the value of the potential coefficient at this mesh
                % index
                vk(vkInd < 1e-5) = potential.Vi(n);
            end
            
            % Create a real space mesh to calculate the potential on
            switch lattice.Dimension
                case '1D'
                    density = 1000;
                    cells = [0 0];
                case'2D'
                    density = [100 100];
                    cells = [0 0; 0 0];
                case '3D'
                    density = [20 20 20];
                    cells = [0 0; 0 0; 0 0];
            end
            superCell = SuperCell(lattice.Dimension, lattice.R, density, cells);
            % The values of the potential on the real space mesh
            v = vk * exp(1i * transpose(recip.Gk) * superCell.Mesh);
            vMax = max(real(v)); vMin = min(real(v));
            % Renormalise the coefficients so that the lattice depth is V0
            vk(kCentre) = vk(kCentre) - vMin;
            vk = potential.V0 * sqrt(lattice.Vol) * vk / (vMax - vMin);
            potential.Vk = vk;
        end
        
        function Plot(potential, lattice, recip, superCell, isoSurf)
            % PLOT THE POTENTIAL
            %
            % REQUIRED INPUT ARGUMENTS
            %   lattice --> instance of the Lattice class
            %     recip --> instance of the ReciprocalLattice class
            % superCell --> instance of the SuperCell class
            %
            % OPTIONAL INPUT ARGUMENTS
            % isosurf --> value of the iso-surface for a 3D plot
            
            if nargin < 5
                isoSurf = 0; end
            % Calculate the potential at the real space mesh coordinates using the Fourier components of
            % the potential
            potentialMesh = potential.Vk * exp(1i * transpose(recip.Gk) * superCell.Mesh) / sqrt(lattice.Vol);
            potentialMesh = real(potentialMesh);
            isoSurf = isoSurf * potential.V0;
            % Plot the potential
            figure;
            switch lattice.Dimension
                case '1D'
                    plot(superCell.MeshX, potentialMesh);
                    axis([min(superCell.MeshX) max(superCell.MeshX) 0 potential.V0]);
                    xlabel('x/\lambda'); ylabel('V/E_R');
                    title('Optical lattice potential');
                case '2D'
                    potentialMesh = reshape(potentialMesh, superCell.NumXpts, superCell.NumYpts);
                    contours = linspace(0, potential.V0, 11);
                    contourf(superCell.MeshX, superCell.MeshY, potentialMesh, contours);
                    axis([min(superCell.MeshX) max(superCell.MeshX) min(superCell.MeshY) max(superCell.MeshY)]);
                    xlabel('x/\lambda'); ylabel('y/\lambda'); zlabel('V/E_R');
                    title('Optical lattice potential');
                    colorbar;
                    caxis([0,potential.V0]);
                    
                case '3D'
                    potentialMesh = reshape(potentialMesh, superCell.NumXpts, superCell.NumYpts, ...
                        superCell.NumZpts);
                    fv = patch(isosurface(superCell.MeshX, superCell.MeshY, superCell.MeshZ, ...
                        potentialMesh, isoSurf));
                    axis([min(superCell.MeshX) max(superCell.MeshX) min(superCell.MeshY) ...
                        max(superCell.MeshY) min(superCell.MeshZ) max(superCell.MeshZ)]);
                    xlabel('x/\lambda'); ylabel('y/\lambda'); zlabel('z/\lambda');
                    title('Optical lattice potential');
                    set(fv, 'FaceColor', 'blue', 'EdgeColor', 'none')
                    view([-65, 20])
                    axis tight
                    camlight left;
                    set(gcf, 'Renderer', 'zbuffer'); lighting phong
            end
        end
    end
end
