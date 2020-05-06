function [] = Wannier(G,hkl,vi,v0)

if (nargin < 1)
        G = 4 * pi * [[1; -1] [1; 1]];
        vi = [0.2395+0.0342i 0.2395-0.0342i 0.0570+0.0434i 0.0570-0.0434i 0.2001i -0.2001i 0.0570-0.1567i 0.0570+0.1567i];
        hkl = [[1; -1] [-1; 1] [1; 1] [-1;-1] [2; 0] [-2; 0] [0; 2] [0; -2]];
        v0 = 30;
        vi = -vi;
        %try to do the 30 Er to see localization in each well.
end
% Performs a single run of RunWannier.m and plots the potential and Wannier functions
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Giovanni Cotugno, Tomi Johnson, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

% #########################################################################################
% INPUT PARAMETERS
% #########################################################################################

% ---------------------
% THE POTENTIAL
% ---------------------

% Please decide on some reference length scale lambda. All units of inputs
% and outputs will be in terms of this length.
% The length scale could, for example, be the wavelength of the lasers
% creating your optical lattice or the lattice parameter of the lattice you
% wish to create.

% Primitive reciprocal lattice vectors. Lengths are given in units of 1/lambda
% G = (1./sqrt(2))* 4 * pi * [[1; -1] [1; 1]];
% G = 4*pi*[[1; 0] [0; 1]];
% The coordinates of reciprocal lattice vectors corresponding to non-zero 
% potential coefficients, given in units of the primitive reciprocal 
% lattice vectors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hkl = [[1; -1] [1; 1] [-1; 1] [-1;-1] [2; 0] [0; 2] ...
%        [-1; 1] [-1; -1] [1; -1] [1; 1] [-2; 0] [0; -2]];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hkl = [[1; -1] [-1; 1] [1; 1] [-1;-1] [2; 0] [-2; 0] [0; 2] [0; -2]];
% hkl = [[0; 0] [1; 0] [-1; 0] [0; 1] [0; -1]];
% vi = [1 -1 -1 -1 -1];
% The corresponding potential coefficients up to some positive constant of
% proportionality (their magnitude will be determined by the next parameter v0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vi = 0.5.*[2 1.2*exp(1i*90*(pi./180)) 0.6*exp(-1i*160*(pi./180)) 1*exp(1i*70*(pi./180)) 1.2*exp(1i*90*(pi./180)) 1*exp(-1i*70*(pi./180))...
%     2 1.2*exp(-1i*90*(pi./180)) 0.6*exp(1i*160*(pi./180)) 1*exp(-1i*70*(pi./180)) 1.2*exp(-1i*90*(pi./180)) 1*exp(1i*70*(pi./180))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vi = [0.2395+0.0342i 0.2395-0.0342i 0.0570+0.0434i 0.0570-0.0434i 0.2001i -0.2001i 0.0570-0.1567i 0.0570+0.1567i];
% Depth of the potential in real space, i.e. V_max - V_min, in units of the
% recoil energy, E_R = h^2 / 2 mu lambda^2, where mu is the mass of
% the particles
% v0 = 10;

% ---------------------
% APPROXIMATIONS MADE IN THE CALCULATION OF THE WANNIER STATES
% ---------------------

% Which bands would you like to use to calculate the Wannier functions? 
% Also, which bands would you like to be be allowed to mix with each other?
% Example A: groups = {1,2,3}; Calculate the Wannier functions for the lowest three bands. 
% No interband mixing is allowed, corresponding to the calculation of the 
% maximally localized ordinary Wannier functions.
% Example B: groups = {[1,2,3]}; Calculate the Wannier functions for the lowest three bands. 
% Mixing is allowed between all bands, corresponding to the calculation of 
% the maximally localized generalized Wannier functions.
% Example C: groups = {1,[2,3]}; Calculate the Wannier functions for the lowest three bands. 
% Mixing is allowed between bands 2 and 3 only, corresponding to the calculation of 
% the maximally localized ordinary Wannier function for the lowest band,
% and the maximally localized generalized Wannier function for the two excited bands.

groups = {[1,2]};

% The maximum norm of reciprocal lattice vectors included in the mesh for the Fourier
% decompositions. In units of 1/lambda.
Gmax = 100;

% No. of real space unit cells / no. of mesh points in each direction of the 1st BZ.
% The value must be even
N = 50;

% ---------------------
% HUBBARD PARAMETERS OF INTEREST
% ---------------------

% The lattice sites at which to calculate Wannier functions, given in units of the lattice
% vectors. Hopping and interaction parameters will be calculated between Wannier functions
% at each site and at the origin
mbSites = [[0; 0] [1; 0] [0; 1] [1; 1] [1; -1]];

% ---------------------
% APPROXIMATIONS MADE IN THE CALCULATION OF THE HUBBARD INTERACTION
% PARAMETERS
% ---------------------

% mbDensity: density of real space points (per cell) in each lattice direction
mbDensity = [10 10];
% mbCells: no of cells in each direction to produce the real space grid
mbCells = [-5 5; -5 5];

% --------------
% PLOT PARAMETERS
% --------------
% plotCells: limits of the super-cell on which to plot the functions
% e.g. [-1 1] produces a plot including the central unit cells and one
% either side
% plotDensity: no. of points in each real space unit cell, in each lattice direction

plotCells = [-1 1; -1 1];
plotDensity = [25 25];

% ---------------
% CALCULATION ADMIN
% ---------------
filename = 'UniqueName'; % A filename for saving the data (make this parameter dependent for data runs)

recalcBloch = 'true'; % Set true to recalculate the Bloch data
reloadComposite = 'false'; % Set true to load previous wannier data
disentangle = 'true'; % Set true to disentangle the bands before running Wannier90
recalcComposite = 'true'; % Set true to recalculate the wannier90 optimised Bloch data
recalcManyBody = 'true'; % Set true to recalculate the Many Body data

parallelTransport = 'false'; % 1D only. Set true to use the Parallel Transport algorithm
calcMode = 'multiprod'; % 'multiprod' is memory intensive. Use 'loop' if this causes issues 
randomise = 'true'; % Set true to randomise the band indices at each k before disentangling

iterIso = 1000; % The number of iterations for optimising isolated bands
iterComp = 1000; % The number of iterations for optimising composite bands
iterDis = 1000; % The number of iterations for disentangling the bands

epsilon = 0.1; % Variable used in Wannier90. Set between 0 and 1
alphaDis = 1; % Variable used in the disentangling algorithm. Set between 0 and 1

% #########################################################################################
% END OF INPUT PARAMETERS
% #########################################################################################






















% #########################################################################################
% WANNIER CALCULATION (Stages 1-4)
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
% PLOTTING (Stage 5)
% #########################################################################################

% --------------------------------------------------------
% Produce the plots of the potential and Wannier functions
% --------------------------------------------------------
Plots(lattice, potential, recip, bloch, wannier90, manyBody, plotDensity, plotCells, calcMode)
end