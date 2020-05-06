
% A class defining the real space and reciprocal lattice vectors of a
% lattice given a dimension, bravais structure, and lattice parameters, or
% alternatively by specifying the reciprocal lattice vectors
% 
% PROPERTIES
% Dimension: The dimensions of the lattice. '1D', '2D', or '3D'.
% Bravais: Type of Bravais lattice
% A,B,C: Lattice constants for the conventional unit cell
% Chi: Angle between conventional lattice vectors R_1 and R_2
% Theta: Angle between conventional lattice vectors R_1 and R_3
% Phi: Angle between conventional lattice vectors R_2 and R_3
% R: Primitive lattice vectors
% G: Reciprocal lattice vectors
% Vol: Volume of the primitive unit cell
% 
% METHODS:
% Lattice(G): Constructor
% SetVectorsFromG(G): Sets the lattice vectors given the reciprocal lattice vectors
% SetVectorsFromBravais(dimension, bravais, a, b, c, chi, theta, phi): Sets
% the lattice and reciprocal lattice vectors given a Bravais structure and parameters
% MinimalGSet(lattice): Finds the set of reciprocal lattice vectors with the smallest norms
% 
% HIDDEN METHODS
% CheckInputs(dimension, bravais): Checks whether dimension and bravais are of the correct form
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

classdef Lattice
    
    properties(GetAccess = 'public', SetAccess = 'private')
        % DECLARE PROPERTIES AND DEFAULT VALUES
        Dimension;
        Bravais = 'line';
        A;
        B = 0;
        C = 0;
        Chi = 90;
        Theta = 90;
        Phi = 90;
        R;
        G;
        Vol;
    end
    
    methods
        function lattice = Lattice(G)
            % CLASS CONSTRUCTOR
            %
            % REQUIRED INPUT ARGUMENTS
            %       G --> reciprocal lattice vectors, e.g. [1, 0; 0, 1]
            
            if nargin == 1
                lattice = lattice.SetVectorsFromG(G); end
        end
        
        function lattice = SetVectorsFromG(lattice, G)
            
            lattice.G = G;
            switch size(G, 1)
                case 1
                    lattice.Dimension = '1D';
                case 2
                    lattice.Dimension = '2D';
                case 3
                    lattice.Dimension = '3D';
            end
            
            % Set the volume and lattice vectors
            switch lattice.Dimension
                case '1D'
                    lattice.R = 2 * pi / lattice.G;
                    lattice.Vol = lattice.R;
                case '2D'
                    lattice.G(3,:) = [0 0]; lattice.G(:,3) = [0; 0; 1];
                    recipVol = det(lattice.G');
                    lattice.R(:,1) = (2 * pi / recipVol) * cross(lattice.G(:,2), lattice.G(:,3));
                    lattice.R(:,2) = (2 * pi / recipVol) * cross(lattice.G(:,3), lattice.G(:,1));
                    lattice.R(:,3) = (2 * pi / recipVol) * cross(lattice.G(:,1), lattice.G(:,2));
                    lattice.R = lattice.R(1:2,1:2); lattice.G = lattice.G(1:2,1:2);
                    lattice.Vol = det(lattice.R');
                case '3D'
                    recipVol = det(lattice.G');
                    lattice.R(:,1) = (2 * pi / recipVol) * cross(lattice.G(:,2), lattice.G(:,3));
                    lattice.R(:,2) = (2 * pi / recipVol) * cross(lattice.G(:,3), lattice.G(:,1));
                    lattice.R(:,3) = (2 * pi / recipVol) * cross(lattice.G(:,1), lattice.G(:,2));
                    lattice.Vol = det(lattice.R');
            end
        end
        
        function lattice = SetVectorsFromBravais(lattice, dimension, bravais, a, b, c, chi, theta, phi)
            % CLASS CONSTRUCTOR
            %
            % REQUIRED INPUT ARGUMENTS
            %       dimension --> lattice.Dimension    '1D', '2D', or '3D'
            %         bravais --> lattice.Bravais      e.g. 'square'
            %           a,b,c --> lattice.A,B,C        e.g. 1
            % chi, theta, phi --> lattice.Cells        e.g. 90
            
            Lattice.CheckInputs(dimension, bravais);
            lattice.Dimension = dimension;
            if strcmp(dimension, '2D') || strcmp(dimension, '3D')
                lattice.Bravais = bravais; end
            lattice.A = a;
            
            % Set the lattice parameters
            switch dimension
                case '1D'
                    lattice.R = a;
                case '2D'
                    switch bravais
                        case 'square'
                            lattice.B = a;
                            lattice.R(:,1) = [a; 0];
                            lattice.R(:,2) = [0; a];
                        case 'rectangular-P'
                            lattice.B = b;
                            lattice.R(:,1) = [a; 0];
                            lattice.R(:,2) = [0; b];
                        case 'rectangular-C'
                            lattice.B = b;
                            lattice.R(:,1) = [a; 0];
                            lattice.R(:,2) = [a/2; b/2];
                        case 'hexagonal'
                            lattice.B = a;
                            lattice.Chi = 60;
                            lattice.R(:,1) = [a; 0];
                            lattice.R(:,2) = [a/2; a*sind(60)];
                        case 'oblique'
                            lattice.B = b;
                            lattice.Chi = chi;
                            lattice.R(:,1) = [a; 0];
                            lattice.R(:,2) = [b*cosd(chi); b*sind(chi)];
                    end
                case '3D'
                    switch bravais
                        case 'cubic-P'
                            lattice.B = a; lattice.C = a;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [0; a; 0];
                            lattice.R(:,3) = [0; 0; a];
                        case 'cubic-I'
                            lattice.B = a; lattice.C = a;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [0; a; 0];
                            lattice.R(:,3) = [a/2; a/2; a/2];
                        case 'cubic-F'
                            lattice.B = a; lattice.C = a;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [a/2; a/2; 0];
                            lattice.R(:,3) = [a/2; 0; a/2];
                        case 'tetragonal-P'
                            lattice.B = b; lattice.C = a;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [0; b; 0];
                            lattice.R(:,3) = [0; 0; a];
                        case 'tetragonal-I'
                            lattice.B = b; lattice.C = a;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [0; b; 0];
                            lattice.R(:,3) = [a/2; b/2; a/2];
                        case 'orthorhombic-P'
                            lattice.B = b; lattice.C = c;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [0; b; 0];
                            lattice.R(:,3) = [0; 0; c];
                        case 'orthorhombic-C'
                            lattice.B = b; lattice.C = c;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [a/2; b/2; 0];
                            lattice.R(:,3) = [0; 0; c];
                        case 'orthorhombic-I'
                            lattice.B = b; lattice.C = c;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [0; b; 0];
                            lattice.R(:,3) = [a/2; b/2; c/2];
                        case 'orthorhombic-F'
                            lattice.B = b; lattice.C = c;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [a/2; b/2; 0];
                            lattice.R(:,3) = [a/2; 0; c/2];
                        case 'rhombohedral'
                            lattice.B = a; lattice.C = a;
                            lattice.Chi = chi; lattice.Theta = chi; lattice.Phi = chi;
                            cosphi2 = (cosd(chi) - cosd(chi) * cosd(chi)) / (sind(chi) * sind(chi));
                            sinphi2 = sqrt(sind(chi)^2 - cosd(chi)^2 - cosd(chi)^2 + ...
                                2 * cosd(chi) * cosd(chi) * cosd(chi)) / (sind(chi) * sind(chi));
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = a * [abs(cosd(chi)); abs(sind(chi)); 0];
                            lattice.R(:,3) = a * [cosd(chi); cosphi2 * sind(chi); sinphi2 * sind(chi)];
                        case 'monoclinic-P'
                            lattice.B = b; lattice.C = c;
                            lattice.Theta = theta;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [0; b; 0];
                            lattice.R(:,3) = c * [cosd(theta); 0; sind(theta)];
                        case 'monoclinic-C'
                            lattice.B = b; lattice.C = c;
                            lattice.Theta = theta;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = [a/2; b/2; 0];
                            lattice.R(:,3) = c * [cosd(theta); 0; sind(theta)];
                        case 'triclinic'
                            lattice.B = b; lattice.C = c;
                            lattice.Chi = chi; lattice.Theta = theta; lattice.Phi = phi;
                            cosphi2 = (cosd(phi) - cosd(theta) * cosd(chi)) / (sind(theta) * sind(chi));
                            sinphi2 = sqrt(sind(chi)^2 - cosd(theta)^2 - cosd(phi)^2 + ...
                                2 * cosd(phi) * cosd(theta) * cosd(chi)) / (sind(theta) * sind(chi));
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = b * [abs(cosd(chi)); abs(sind(chi)); 0];
                            lattice.R(:,3) = c * [cosd(theta); cosphi2 * sind(theta); sinphi2 * sind(theta)];
                        case 'hexagonal'
                            lattice.B = a; lattice.C = c;
                            lattice.Chi = 60;
                            lattice.R(:,1) = [a; 0; 0];
                            lattice.R(:,2) = a * [abs(cosd(60)); abs(sind(60)); 0];
                            lattice.R(:,3) = [0; 0; c];
                    end
            end
            
            % Set the volume and reciprocal lattice vectors
            switch dimension
                case '1D'
                    lattice.Vol = lattice.R;
                    lattice.G = 2 * pi / lattice.R;
                case '2D'
                    lattice.R(3,:) = [0 0]; lattice.R(:,3) = [0; 0; 1]; lattice.Vol = det(lattice.R');
                    lattice.G(:,1) = (2 * pi / lattice.Vol) * cross(lattice.R(:,2), lattice.R(:,3));
                    lattice.G(:,2) = (2 * pi / lattice.Vol) * cross(lattice.R(:,3), lattice.R(:,1));
                    lattice.G(:,3) = (2 * pi / lattice.Vol) * cross(lattice.R(:,1), lattice.R(:,2));
                    lattice.R = lattice.R(1:2,1:2); lattice.G = lattice.G(1:2,1:2);
                case '3D'
                    lattice.Vol = det(lattice.R');
                    lattice.G(:,1) = (2 * pi / lattice.Vol) * cross(lattice.R(:,2), lattice.R(:,3));
                    lattice.G(:,2) = (2 * pi / lattice.Vol) * cross(lattice.R(:,3), lattice.R(:,1));
                    lattice.G(:,3) = (2 * pi / lattice.Vol) * cross(lattice.R(:,1), lattice.R(:,2));
            end
        end
        
        function lattice = MinimalGSet(lattice)
            % FIND THE SET OF RECIPROCAL LATTICE VECTORS WITH THE SMALLEST NORMS
            
            if strcmp(lattice.Dimension, '2D') || strcmp(lattice.Dimension, '3D')
                % Create a reciprocal lattice mesh using current reciprocal lattice basis vectors
                switch lattice.Dimension
                    case '2D'
                        d = [51; 51];
                        numKpts = prod(d); % Total no. of reciprocal lattice points included in DFT
                        [gnx, gny] = meshgrid((d(1) - 1) / 2 : -1 : -(d(1) - 1) / 2, ...
                            (d(2) - 1) / 2 : -1 : -(d(2) - 1) / 2);
                        % Fractional coordinates of the reciprocal lattice vectors
                        gFrac(1,:) = reshape(gnx, 1, numKpts);
                        gFrac(2,:) = reshape(gny, 1, numKpts);
                    case '3D'
                        d = [21; 21; 21];
                        numKpts = prod(d); % Total no. of reciprocal lattice points included in DFT
                        [gnx, gny, gnz] = meshgrid((d(1) - 1) / 2 : -1 : -(d(1) - 1) / 2, ...
                            (d(2) - 1) / 2 : -1 : -(d(2) - 1) / 2, (d(3) - 1) / 2 : -1 : - (d(3) - 1) / 2);
                        % Fractional coordinates of the reciprocal lattice vectors
                        gFrac(1,:) = reshape(gnx, 1, numKpts);
                        gFrac(2,:) = reshape(gny, 1, numKpts);
                        gFrac(3,:) = reshape(gnz, 1, numKpts);
                end
                % The reciprocal lattice vectors included in the DFT, labelled by k
                gk = lattice.G * gFrac;
                % Sort mesh by length
                gNorm = sqrt(sum(gk .^ 2, 1)); [gNorm, gInd] = sort(gNorm);
                % Set the basis vector G1
                lattice.G(:, 1) = gk(:, gInd(2));
                % Find the next smallest basis vector that is linearly independent of G1
                depTest = 1; loop = 3;
                while abs(depTest - 1) < 0.000001 % Stop if dot(G1, G2) / |G1||G2| ~= 1
                    lattice.G(:, 2) = gk(:, gInd(loop)); % New basis vector, G2
                    % Linear dependence test (dot(G1, G2) / |G1||G2| ~ =1)
                    depTest = abs(dot(lattice.G(:,1) / norm(lattice.G(:,1)), lattice.G(:,2) / ...
                        norm(lattice.G(:,2))));
                    loop = loop + 1;
                end
                if strcmp(lattice.Dimension, '3D')
                    % Find the next smallest basis vector that is linearly independent of G1 and G2
                    depTest = 0;
                    while depTest == 0 % Stop if det(G) ~= 0
                        lattice.G(:,3) = gk(:, gInd(loop)); % New basis vector, G2
                        % Linear dependence test (det(G) ~= 0)
                        depTest = det(lattice.G); loop = loop + 1;
                    end
                end
            end
        end
    end
    
    methods(Static = true, Hidden = true, Access = private)
        function CheckInputs(dimension, bravais)
            % CHECK IF THE dimension AND bravais INPUTS ARE OF THE CORRECT FORM
            
            if ~(strcmp(dimension, '1D') || strcmp(dimension, '2D') || strcmp(dimension, '3D'))
                error('Dimension must be set as 1D, 2D, or 3D.'); end
            switch dimension
                % In 2D there are 5 fundamental Bravais lattices
                case '2D'
                    if ~(strcmp(bravais, 'square') || strcmp(bravais, 'rectangular-P') ...
                            || strcmp(bravais, 'rectangular-C') || strcmp(bravais, 'hexagonal') ...
                            || strcmp(bravais, 'oblique'))
                        error(['In 2D, Bravais must be set as square, rectangular-P,' ...
                            ' rectangular-C, hexagonal or oblique.'])
                    end
                % In 3D there are 14 fundamental Bravais lattices
                case '3D'
                    if ~(strcmp(bravais, 'cubic-P') || strcmp(bravais, 'cubic-I') ...
                            || strcmp(bravais, 'cubic-F') ||strcmp(bravais, 'tetragonal-P') ...
                            || strcmp(bravais, 'tetragonal-I') || strcmp(bravais, 'orthorhombic-P') ...
                            || strcmp(bravais, 'orthorhombic-C') || strcmp(bravais, 'orthorhombic-I') ...
                            || strcmp(bravais, 'orthorhombic-F') || strcmp(bravais, 'rhombohedral') ...
                            || strcmp(bravais, 'monoclinic-P') || strcmp(bravais, 'monoclinic-C') ...
                            || strcmp(bravais, 'triclinic') || strcmp(bravais, 'hexagonal'))
                        error(['In 3D, Bravais must be set as cubic-P, cubic-I, cubic-F,' ...
                            ' tetragonal-P, tetragonal-I, orthorhombic-P, orthorhombic-C,' ...
                            ' orthorhombic-I, orthorhombic-F, rhombohedral, monoclinic-P,' ...
                            ' monoclinic-C, triclinic, or hexagonal.'])
                    end
            end
        end
    end
end