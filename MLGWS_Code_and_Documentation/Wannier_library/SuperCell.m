
% Defines a superCell containing mesh points on which to numerically compute functions
%
% PROPERTIES
% Dimension: Dimensions of the superCell 
% R: Primitive lattice vectors for a single cell
% Cells: Upper and lower bounds on the cells in each direction
% NumCells: Number of cells in each direction
% Density: Number of mesh points in each direction per unit cell
% Overlap: Set true to include overlapping points on the cell boundaries
% Mesh: A uniform mesh of points through the superCell.
% MeshX, MeshY, MeshZ: The mesh points in the each direction (useful for plotting functions across the superCell)
% MeshInd: Matrix of indices for the mesh points
% NumPts: Total number of mesh points
% NumXpts, NumYpts, NumZpts: Number of mesh points in each direction
% GammaShift: The shift from the origin of the central mesh point
% dVol: The volume associated with a mesh point
%
% METHODS
% SuperCell(dimension, r, density, cells, gammaShift, overlap): Constructor
% Plot(): Plot the supercell mesh
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

classdef SuperCell
    
    properties(GetAccess = 'public', SetAccess = 'private')
        % DECLARE PROPERTIES AND DEFAULT VALUES
        Dimension;
        R;
        Cells;
        NumCells;
        Density;
        Overlap;
        Mesh;
        MeshX;
        MeshY;
        MeshZ;
        MeshInd;
        NumPts;
        NumXpts;
        NumYpts;
        NumZpts;
        GammaShift ;
        dVol;
    end
    
    methods
        function superCell = SuperCell(dimension, r, density, cells, overlap, gammaShift)
            % CLASS CONSTRUCTOR
            %
            % REQUIRED INPUT ARGUMENTS
            %  dimension --> superCell.Dimension    '1D', '2D', or '3D'
            %          r --> superCell.R            e.g. [1 0; 0 1]
            %    density --> superCell.Density      e.g. [20 20]
            %      cells --> superCell.Cells        e.g. [0 0; -1 1]
            %
            % OPTIONAL INPUT ARGUMENTS
            %    overlap --> superCell.Overlap      default: 'false'
            % gammaShift --> superCell.GammaShift   default: [0 0 0]
            
            % Set properties
            if nargin < 6
                superCell.GammaShift = [0 0 0];
            else
                superCell.GammaShift = gammaShift;
            end
            if nargin < 5
                superCell.Overlap = 'false';
            else
                superCell.Overlap = overlap;
            end
            superCell.Dimension = dimension;
            superCell.R = r;
            switch superCell.Dimension
                case '1D'
                    cells = cells(1,:);
                    density = density(1);
                    superCell.GammaShift = superCell.GammaShift(1);
                case '2D'
                    cells = cells(1:2,:);
                    density = density(1:2);
                    superCell.GammaShift = superCell.GammaShift(1:2);
            end
            superCell.Cells = cells;
            superCell.NumCells = cells(:,2) - cells(:,1) + ones(size(r, 1), 1);
            superCell.Density = density;
            % Set dVol
            switch superCell.Dimension
                case '1D'
                    superCell.dVol = r / density;
                case '2D'
                    r(3,:) = [0 0]; r(:,3) = [0; 0; 1];
                    superCell.dVol = det(r') / prod(density);
                    r = r(1:2,1:2);
                case '3D'
                    superCell.dVol = det(r') / prod(density);
            end
            % Fractional mesh coordinates in each direction
            switch superCell.Dimension
                case '1D'
                    fracMesh = cells(1) - 0.5 : 1 / density(1) : cells(2) + 0.5;
                case '2D'
                    fracMeshX = cells(1,1) - 0.5 : 1 / density(1) : cells(1,2) + 0.5;
                    fracMeshY = cells(2,1) - 0.5 : 1 / density(2) : cells(2,2) + 0.5;
                case '3D'
                    fracMeshX = cells(1,1) - 0.5 : 1 / density(1) : cells(1,2) + 0.5;
                    fracMeshY = cells(2,1) - 0.5 : 1 / density(2) : cells(2,2) + 0.5;
                    fracMeshZ = cells(3,1) - 0.5 : 1 / density(3) : cells(3,2) + 0.5;
            end
            % Remove overlapping mesh points if Overlap != 'true'
            if ~strcmp(superCell.Overlap, 'true')
                switch dimension
                    case '1D'
                        fracMesh = fracMesh(1 : end - 1);
                    case '2D'
                        fracMeshX = fracMeshX(1 : end - 1);
                        fracMeshY = fracMeshY(1 : end - 1);
                    case '3D'
                        fracMeshX = fracMeshX(1 : end - 1);
                        fracMeshY = fracMeshY(1 : end - 1);
                        fracMeshZ = fracMeshZ(1 : end - 1);
                end
            end
            % Set Mesh points in each direction
            switch superCell.Dimension
                case '1D'
                    superCell.MeshX = fracMesh * r;
                case '2D'
                    superCell.MeshX = fracMeshX * norm(r(:, 1));
                    superCell.MeshY = fracMeshY * norm(r(:, 2));
                case '3D'
                    superCell.MeshX = fracMeshX * norm(r(:, 1));
                    superCell.MeshY = fracMeshY * norm(r(:, 2));
                    superCell.MeshZ = fracMeshZ * norm(r(:, 3));
            end
            superCell.NumXpts = length(superCell.MeshX);
            superCell.NumYpts = length(superCell.MeshY);
            superCell.NumZpts = length(superCell.MeshZ);
            % Set the mesh indices
            switch superCell.Dimension
                case '1D'
                    fracMesh = fracMesh + superCell.GammaShift(1);
                    numPts = length(fracMesh);
                    superCell.MeshInd = 1 : numPts;
                case '2D'
                    numPtsTemp = [length(fracMeshX); length(fracMeshY)];
                    numPts = numPtsTemp(1) * numPtsTemp(2);
                    [meshTempX, meshTempY] = meshgrid(fracMeshX, fracMeshY);
                    fracMesh(1,:) = reshape(meshTempX, 1, numPts) + superCell.GammaShift(1);
                    fracMesh(2,:) = reshape(meshTempY, 1, numPts) + superCell.GammaShift(2);
                    superCell.MeshInd = reshape(1 : numPts, numPtsTemp(2), ...
                        numPtsTemp(1));
                case '3D'
                    numPtsTemp = [length(fracMeshX); length(fracMeshY); length(fracMeshZ)];
                    numPts = numPtsTemp(1) * numPtsTemp(2) * numPtsTemp(3);
                    [meshTempX, meshTempY, meshTempZ] = meshgrid(fracMeshX, fracMeshY, fracMeshZ);
                    fracMesh(1,:) = reshape(meshTempX, 1, numPts) + superCell.GammaShift(1);
                    fracMesh(2,:) = reshape(meshTempY, 1, numPts) + superCell.GammaShift(2);
                    fracMesh(3,:) = reshape(meshTempZ, 1, numPts) + superCell.GammaShift(3);
                    superCell.MeshInd = reshape(1 : numPts, numPtsTemp(2), ...
                        numPtsTemp(1), numPtsTemp(3));
            end
            % Set the mesh points
            superCell.Mesh = r * fracMesh;
            superCell.NumPts = numPts;
        end
        
        function Plot(superCell)
            % PLOT THE SUPERCELL MESH
            switch superCell.Dimension
                case '1D'
                    y = zeros(1, superCell.NumPts);
                    plot(superCell.Mesh, y, 'b.');
                case '2D'
                    plot(superCell.Mesh(1,:), superCell.Mesh(2,:), 'b.');
                    set(gca, 'DataAspectRatio', [1 1 1]);
                case '3D'
                    plot3(superCell.Mesh(1,:), superCell.Mesh(2,:), superCell.Mesh(3,:), 'b.');
                    set(gca, 'DataAspectRatio', [1 1 1]);
            end
        end
    end
end