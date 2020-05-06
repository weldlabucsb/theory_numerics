
% Defines a reciprocal lattice mesh given a lattice and an energy cut-off.
% 
% PROPERTIES
% Dimension: The Dimension of the reciprocal lattice
% Size: Initial no. of reciprocal lattice points in each direction
% NumKpts: The total number of reciprocal lattice points
% Gk: An ordered list of reciprocal lattice vectors
% Gfrac: An ordered list of reciprocal lattice vectors given by their fractional coordinate
% Gmax: The maximum norm of reciprocal lattice vectors included in the mesh
% Gkm: Matrix of Fourier coefficient labels
% GkiIn: Index of the original reciprocal lattice points included within the cut-off radius
% 
% METHODS
% ReciprocalLattice(dimension, G, Gmax): Constructor
% DifferenceMatrix(): Constructs a matrix of coefficient labels for difference vectors between
% reciprocal lattice points, with (NumKpts - 1) / 2 as the origin. (Used to calculate v_{k-m}).
% CutOff(): Removes reciprocal lattice points with a norm greater than Gmax
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

classdef ReciprocalLattice
    
    properties(GetAccess = 'public', SetAccess = 'private')
        % DECLARE PROPERTIES AND DEFAULT VALUES
        Dimension;
        Size;
        NumKpts;
        Gk;
        Gfrac;
        Gmax;
        Gkm;
        GkiIn;
    end
    
    methods
        
        function recip = ReciprocalLattice(dimension, G, Gmax)
            % CLASS CONSTRUCTOR
            %
            % REQUIRED INPUT ARGUMENTS
            %    dimension --> recip.Dimension
            %            G --> the primitive reciprocal lattice vectors
            % Gmax --> recip.Gmax
            
            recip.Dimension = dimension;
            recip.Gmax = Gmax;
            switch recip.Dimension
                case '1D'
                    recip.Size = 2 * floor(Gmax / G) + 1;
                    recip.NumKpts = recip.Size;
                    recip.Gfrac = -(recip.Size - 1) / 2 : (recip.Size - 1) / 2;
                case '2D'
                    gTemp = diag(ones(1,3)); gTemp(1:2,1:2) = G; % Add the z-direction
                    recip.Size(1) = 2 * floor(abs(Gmax * norm(cross(gTemp(:,2), gTemp(:,3)))/ ...
                        dot(gTemp(:,1), cross(gTemp(:,2), gTemp(:,3))))) + 1;
                    recip.Size(2) = 2 * floor(abs(Gmax * norm(cross(gTemp(:,1), gTemp(:,3)))/ ...
                        dot(gTemp(:,2), cross(gTemp(:,1), gTemp(:,3))))) + 1;
                    recip.NumKpts = prod(recip.Size);
                    [gnx, gny] = meshgrid(-(recip.Size(1) - 1) / 2 : (recip.Size(1) - 1) / 2, ...
                        -(recip.Size(2)-1) / 2 : (recip.Size(2)-1) / 2);
                    recip.Gfrac(1,:) = reshape(gnx, 1, recip.NumKpts);
                    recip.Gfrac(2,:) = reshape(gny, 1, recip.NumKpts);
                case '3D'
                    recip.Size(1) = 2 * floor(abs(Gmax * norm(cross(G(:,2), G(:,3))) / ...
                        dot(G(:,1), cross(G(:,2),G(:,3))))) + 1;
                    recip.Size(2) = 2 * floor(abs(Gmax * norm(cross(G(:,1), G(:,3))) / ...
                        dot(G(:,2), cross(G(:,1), G(:,3))))) + 1;
                    recip.Size(3) = 2 * floor(abs(Gmax * norm(cross(G(:,1), G(:,2))) / ...
                        dot(G(:,3), cross(G(:,1), G(:,2))))) + 1;
                    recip.NumKpts = prod(recip.Size);
                    [gnx, gny, gnz] = meshgrid(-(recip.Size(1) - 1) / 2 : (recip.Size(1) - 1) / 2, ...
                        -(recip.Size(2) - 1) / 2 : (recip.Size(2) - 1) / 2, ...
                        -(recip.Size(3) - 1) / 2 : (recip.Size(3) - 1) / 2);
                    recip.Gfrac(1,:) = reshape(gnx, 1, recip.NumKpts);
                    recip.Gfrac(2,:) = reshape(gny, 1, recip.NumKpts);
                    recip.Gfrac(3,:) = reshape(gnz, 1, recip.NumKpts);
            end
            recip.Gk = G * recip.Gfrac;
        end
        
        function recip = DifferenceMatrix(recip)
            % Constructs a matrix of coefficient labels for difference vectors between
            % reciprocal lattice points, with (NumKpts - 1) / 2 as the origin. (Used to calculate v_{k-m})
            
            gkm = diag(ones(recip.NumKpts, 1) * ((recip.NumKpts + 1) / 2));
            for loop = 1:(recip.NumKpts - 1) / 2
                gkm = gkm + diag(ones(recip.NumKpts - loop, 1) * ((recip.NumKpts + 1) / 2 - loop), loop);
                gkm = gkm + diag(ones(recip.NumKpts - loop, 1) * ((recip.NumKpts + 1) / 2 + loop), -loop);
            end
            if strcmp(recip.Dimension, '2D') || strcmp(recip.Dimension, '3D')
                gkm2 = diag(ones(1, recip.Size(2)));
                for loop = 1:(recip.Size(2) - 1) / 2
                    gkm2 = gkm2 + diag(ones(1, recip.Size(2) - loop), loop) + ...
                        diag(ones(1, recip.Size(2) - loop), -loop);
                end
                gkm2 = repmat(gkm2, recip.Size(1), recip.Size(1));
                if strcmp(recip.Dimension, '3D')
                    gkm3 = diag(ones(1, recip.Size(1) * recip.Size(2)));
                    for loop = 1 : (recip.Size(1) * recip.Size(2) - 1) / 2
                        gkm3 = gkm3 + diag(ones(1, recip.Size(1) * recip.Size(2) - loop), loop) + ...
                            diag(ones(1, recip.Size(1) * recip.Size(2) - loop), -loop);
                    end
                    gkm2 = gkm2 .* gkm3;
                    gkm2 = repmat(gkm2, recip.Size(3), recip.Size(3));
                end
                gkm = gkm .* gkm2;
            end
            recip.Gkm = gkm;
        end
        
        function recip = CutOff(recip)
            % Removes reciprocal lattice points with an energy greater than Gmax
            
            if strcmp(recip.Dimension, '2D') || strcmp(recip.Dimension, '3D')
                gkNorm = zeros(1, recip.NumKpts);
                for k = 1 : recip.NumKpts
                    gkNorm(k) = norm(recip.Gk(:,k)); end
                recip.Gk = recip.Gk(:, gkNorm <= recip.Gmax);
                % Index of reciprocal lattice vectors within E_cut
                recip.GkiIn = 1 : recip.NumKpts;
                recip.GkiIn = recip.GkiIn(:, gkNorm <= recip.Gmax);
                % Index of reciprocal lattice vectors outside E_cut
                gkiOut = 1 : recip.NumKpts; gkiOut = gkiOut(:, gkNorm > recip.Gmax);
                % Remove indices pointing to k-points outside E_cut
                for k = 1 : length(gkiOut)
                    recip.Gkm(recip.Gkm == gkiOut(k)) = 0; end
                recip.NumKpts = length(recip.GkiIn); % Total no. of k-points within E_cut
                % Remove rows and columns for k-points outside E_cut
                recip.Gkm = recip.Gkm(recip.GkiIn, recip.GkiIn);
                % Re-index the k-points in the Gkm matrix
                for k = 1:recip.NumKpts
                    recip.Gkm(recip.Gkm == recip.GkiIn(k)) = k; end
                recip.Gkm = sparse(recip.Gkm);
            end
        end
    end
end