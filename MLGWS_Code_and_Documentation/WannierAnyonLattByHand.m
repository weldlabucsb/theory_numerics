% Performs a single run of RunWannier.m and plots the potential and Wannier functions
% Example 8: 2D Kagome lattice
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Giovanni Cotugno, Tomi Johnson, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------
clear; clc; close all;

% #########################################################################################
% INPUT PARAMETERS
% #########################################################################################


%First, getting lattice depth from code
%THESE ARE THE PARAMETERS USED TO GET WHAT IS BELOW
A = [1,1,0.6,0.5];
ph_deg = [0, 0, 90, -70];
th_deg = [0,90,180,270];
pol_deg = [0,0,0,0];
plots = 1;

%slightly deeper hopping well
% A = [1,1,0.6,0.5];
% ph_deg = [0, 0, 95, -80];
% th_deg = [0,90,180,270];
% pol_deg = [0,0,0,0];
% plots = 1;

[waveAmplitudes,deltaKsUnitless,deltaPhis,maxAmp]=GeneralLatticeComplexRepWithComponents(A,ph_deg,th_deg,pol_deg,plots);



% ---------------------
% MAIN INPUT PARAMETERS
% ---------------------
% The reciprocal lattice vectors. Lengths are given in units of wavelength, i.e. lambda = 1
G = 2 * pi * [[1; -1] [1;1]];
% The coordinates of the potential coefficients in the reciprocal lattice, given in units of
% the reciprocal lattice vectors, G
hkl = [[1;0],[0;1],[1;1],[-1;1],...
    [-1;0],[0;-1],[-1;-1],[1;-1]];
% The ratios of the corresponding potential coefficients
vi1 = -[(1+0.3*exp(1i*(8/9)*pi)),  (0.6*exp(1i*(1/2)*pi)+0.5*exp(-1i*(7/18)*pi)), (0.6*1i), (0.5*exp(-1i*(7/18)*pi))];
vi2 = conj(vi1);
% vi = horzcat(vi1,vi2)/maxAmp;
vi = horzcat(vi1,vi2);  % I think you don't actually need to normalize with maxAmp.  It will do that for you.
% Strength of the potential in units of the recoil energy, E_R
v0 = 10; 
% The composite groups to calculate the generalised Wannier functions for.
% Example: groups = {1,[2,3]}; Band 1 is isolated, while bands 2 and 3 form a composite group
groups = {[1,2],3};  % PLAYED WITH WHAT GROUPINGS  -  Mixing bands 1 and 2 works ok, but sometimes, it doesn't converge
% The lattice sites at which to calculate Wannier functions, given in units of the lattice
% vectors. Hopping and interaction parameters will be calculated between Wannier functions
% at each site and at the origin
mbSites = [[0; 0] [1; 0] [0; 1] [1; 1] [1; -1]];

% ---------------
% FILE PARAMETERS
% ---------------
filename = '2D-Anyon-Lattice'; % A filename for saving the data (make this parameter dependent for data runs)

% ----------------
% BLOCH PARAMETERS
% ----------------
recalcBloch = 'true'; % Set true to recalculate the Bloch data
% The maximum norm of reciprocal lattice vectors included in the mesh for the Fourier
% decompositions.
Gmax = 50; %INCREAED THIS
% No. of real space unit cells / no. of q points in each direction of the 1st BZ.
% The value must be even
N = 20; %INCREAED THIS

% --------------------
% WANNIER90 PARAMETERS
% --------------------
recalcComposite = 'true'; % Set true to recalculate the wannier90 optimised Bloch data (and onwards)
reloadComposite = 'false'; % Set true to load previous wannier data (recalcComposite must be true)
parallelTransport = 'false'; % 1D only. Set true to use the Parallel Transport algorithm
disentangle = 'true'; % Set true to disentangle the bands before running Wannier90
randomise = 'true'; % Set true to randomise the band indices at each q before disentangling
iterIso = 1000; % The number of iterations for optimising isolated bands
iterComp = 1000; % The number of iterations for optimising composite bands
iterDis = 1000; % The number of iterations for disentangling the bands
epsilon = 1; % Variable used in Wannier90. Set between 0 and 1
alphaDis = 1; % Variable used in the disenatngling algorithm. Set between 0 and 1

% ------------------------------
% HUBBARD (MANY BODY) PARAMETERS
% ------------------------------
recalcManyBody = 'true'; % Set true to recalculate the Many Body data
calcMode = 'multiprod'; % 'multiprod' is memory intensive. Use 'loop' if this causes issues 
% mbDensity: density of real space points in each lattice direction
mbDensity = [10 10];
% mbCells: no of cells in each direction to produce the real space grid
mbCells = [-5 5; -5 5];

% --------------
% Plot paramters
% --------------
% plotCells: limits of the super-cell on which to plot the functions
% e.g. [-1 1] produces a plot including the unit cells centred at -lambda, 0, lambda
% plotDensity: no. of points in each real space unit cell, in each lattice direction
plotCells = [-2 2; -2 2];
plotDensity = [20 20];

% #########################################################################################
% WANNIER CALCULATION
% #########################################################################################

% -------------------
% PRE-PROCESSING CODE
% -------------------
% Add library folder to the Matlab path
path('./Wannier_library', path);
% Create data directories if they don't exist
if exist('./Wannier_data', 'dir') ~= 7
    mkdir('./Wannier_data'); end
if exist('./Wannier_data/Parameters', 'dir') ~= 7
    mkdir('./Wannier_data/Parameters'); end
% Calculate total number of bands
[numBands] = TotalBands(groups);
% Save the parameters to a file. This is then loaded by RunWannier.m to calculate the
% Wannier functions
save(['./Wannier_data/Parameters/' filename]);

% -----------------------------------------------------------------------------
% CALCULATE THE BLOCH STATES, THE WANNIER FUNCTIONS, AND THE HUBBARD PARAMETERS
% -----------------------------------------------------------------------------
[lattice, recip, potential, bloch, neighbours, wannier90, manyBody] = RunWannier(filename);

% #########################################################################################
% PLOTTING
% #########################################################################################

% --------------------------------------------------------
% Produce the plots of the potential and Wannier functions
% --------------------------------------------------------
Plots(lattice, potential, recip, bloch, wannier90, manyBody, plotDensity, plotCells, calcMode)
