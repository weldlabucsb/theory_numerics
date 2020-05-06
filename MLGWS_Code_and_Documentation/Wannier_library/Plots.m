
% Plot the potential and the Wannier functions
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

function Plots(lattice, potential, recip, bloch, wannier90, manyBody, plotDensity, plotCells, calcMode)

disp('-------------------------------------------------');
disp('MAXIMALLY LOCALISED GENERALISED WANNIER FUNCTIONS');
disp('Plotting the potential and Wannier functions...');
disp('-------------------------------------------------');
switch lattice.Dimension
    case '1D'
        r = 1;
    case '2D'
        r = [1 0; 0 1];
    case '3D'
        r = [1 0 0; 0 1 0; 0 0 1];
end
% Create a superCell with orthogonal axes
superCell = SuperCell(lattice.Dimension, r, plotDensity, plotCells, 'true');
% Create an instance of ManyBody for the orthogonal superCell
plotWannier = ManyBody(wannier90.Groups, manyBody.Sites, superCell);
% Calculate the Wannier functions
plotWannier = plotWannier.CalculateWannier(lattice, recip, bloch, wannier90, calcMode);
% Plot the potential
potential.Plot(lattice, recip, superCell);
% Plot the Wannier functions
plotWannier.Plot(1, 0.5);