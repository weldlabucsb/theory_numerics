
% Calculate the total number of bands used in the calculations
% 
% -----------------------------------------------------------------------------------------
% Richard Walters, Stephen Clark and Dieter Jaksch.
% Atomic and Laser Physics, Clarendon Laboratory, University of Oxford, Oxford, OX1 3PU, UK
% -----------------------------------------------------------------------------------------

function [numBands] = TotalBands(groups)

numBands = 0;
for n = 1 : length(groups)
    for m = 1 : length(groups{n})
        numBands = numBands + 1;
    end
end
