%% metaWannier.m
%The purpose of this code is to run repeated iterations of Wannier.m and to
%save the results to plot things like the tunnelling element versus lattice
%strength etc...

%% Initialize Lattice As Interfering Beams
%In general the code needs to have everything that the Wannier code would
%have as an input, except in vector form that is iterated over.
%regular lattice
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

deltaKsUnitless = round(deltaKsUnitless);
potentialDepth = 10; %in Er!!
waveAmplitudes = waveAmplitudes.*(potentialDepth./maxAmp);

%% Find the Complex Exponential Coefficients
max_m = 2;
mLength = 2*max_m + 1;
Vcoeff = zeros(mLength,mLength);
for jj = 1:length(waveAmplitudes)
    xKcomp = deltaKsUnitless(jj,1);
    yKcomp = deltaKsUnitless(jj,2);
    %I am making the index xKcomp+(max_m+1)because x/y K components could
    %be negative. Therefore, I want to make sure that the 'center' of the
    %matrix corresponds to the (0,0) coefficient. This is also related to
    %my above comment, since if the basis is too small, then these
    %arguments could be negative and matlab will complain :(
    Vcoeff(xKcomp+(max_m+1),yKcomp+(max_m+1)) = Vcoeff(xKcomp+(max_m+1),yKcomp+(max_m+1)) + waveAmplitudes(jj).*(exp(-1i*deltaPhis(jj)))./2;
    Vcoeff(-xKcomp+(max_m+1),-yKcomp+(max_m+1)) = Vcoeff(-xKcomp+(max_m+1),-yKcomp+(max_m+1)) + waveAmplitudes(jj).*(exp(1i*deltaPhis(jj)))./2;
end
Vcoeff = -Vcoeff;
hkl = [];
vi = [];
for ii = -max_m:max_m
    for jj = -max_m:max_m
        if (Vcoeff(ii + (max_m+1),jj + (max_m+1)) ~= 0)
            hkl = [hkl [ii ; jj]];
            vi = [vi Vcoeff(ii + (max_m+1),jj + (max_m+1))];
        end
    end
end

%% Specify Inputs to the Wannier 90 Code
%We need to put everything in terms of the reciprocal lattice parameters
G = 4 * pi * [[1; -1] [1; 1]];
v0 = 10;
% vi = [0.2395+0.0342i 0.2395-0.0342i 0.0570+0.0434i 0.0570-0.0434i 0.2001i -0.2001i 0.0570-0.1567i 0.0570+0.1567i];
% hkl = [[1; -1] [-1; 1] [1; 1] [-1;-1] [2; 0] [-2; 0] [0; 2] [0; -2]];
%What we want to iterate over:


%% Re-Express in the 45 degree tilted basis
if (0)
    hklp = [];
    G = 4*pi*[[1;0] [0;1]];
    a = [[1 1];[1 -1]];
    for jj = 1:length(hkl)
        hklp = [hklp a\hkl(:,jj)];
    end
    hkl = hklp;
end
keyboard;



% Wannier(G,hkl,vi,v0);
v0 = linspace(10,10,5);
Js = zeros(6,length(v0));
for ii = 1:length(v0)
    Wannier(G,hkl,vi,v0(ii));
    data = load('.\Wannier_data\Many_body_data\UniqueName.mat')
    data = data.manyBody;
    J = data.J;
%     J_01 = J(1,1,2);
%     Js(ii) = J_01;
    Js(1,ii) = J(2,1,1);
    Js(2,ii) = J(2,2,2);
    Js(3,ii) = J(2,2,4);
    Js(4,ii) = J(1,1,2);
    Js(5,ii) = J(1,1,4);
    Js(6,ii) = J(2,1,4);
end
keyboard;
figure;
hold on;
for ii = 1:6
    plot(real(Js(ii,:)))
end
legend('j1','j2','j3','j4','j5','j6');
xlabel('run number, all at V0 = 10 Er');
ylabel('Real(J)');
figure;
hold on;
for ii = 1:6
    plot(imag(Js(ii,:)),v0)
end
legend('j1','j2','j3','j4','j5','j6');
xlabel('run number, all at V0 = 10 Er');
ylabel('Im(J)');
plot(v0,Js);
    