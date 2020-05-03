close all; clear;
%% Constants
% Define physical constants
amu=1.66E-27;       % 1 AMU
m=7*amu;            % Lithium mass
lambda=1064E-9;     % Wavelength of light
h=6.626E-34;        % Planck's Const.
kL=2*pi/lambda;     % Wave Vector
hbar=h/(2*pi);        % Reduced planck's constant
Er=hbar^2*kL^2/(2*m); % Recoil Energy
vR=hbar*kL/m;         % recoil velocity
d=lambda/2;         % Lattice site distance
fr=Er/h;            % recoil frequency in Hz

%% Initialize Parameters
% Drive Parameters
V0 = 4.3;             % average lattice depth Er
alpha = 1/V0;         % percent modulation
fD = 53.56;           % drive frequency Hz
TD = 1/fD;            % drive period s
phi = pi;             % initial drive phase rad
dfdt = 0;             % linear chirp rate

% Force Parameters
w = 2*pi*15.5;           % Harmonic Potential frequency Hz
TB = 16.75E-3;        % Bloch period s
fB = 1/TB;            % Bloch frequency Hz
F = h*fB/d;           % Force Newtons

% BEC Parameters
xwidth = 50*d;        % BEC spatial width in m
kwidth = .5;          % momentum width kL

% Simulation Parameters
xpoints = 11;         % number of position initial conditions (odd)
kpoints = 11;         % number of momentum initial conditions (odd)
xend = 3*xwidth;      % sample width in m
kend = 2*kwidth;      % initial k spacing in kL
tpre = -.01E-3;        % time before modulation starts s
tstep = 1*TD;           % time spacing of saved solutions s
tend = 200*TD;        % simulation end time s
time = 0:tstep:tend;  % time mesh
opts = odeset('Abstol',1E-5,'Reltol',1E-3);

%% Get tunneling/band functions
bw = open('bandwidth.mat');             % contains precomputed bandwidth vs lattice depth values
J = @(V) interp1(bw.depth,bw.BW,V)/4;   % tunneling Er
Ek = @(V,k) -2*J(V)*cos(pi*k);          % ground band Er

%% Set Equations of Motion (1:x in d,2:k in kL, t in 1/fD)
% LATTICE DEPTH over time in Er
Vt = @(t) V0*(1+alpha*sin(2*pi*(1+dfdt*t/fD^2).*t+phi)).*(t>0)+...
    V0*(1+alpha*sin(phi)).*(t<=0);

% GROUP VELOCITY in lattice sites/drive period
vg = @(V,k) J(V)*pi*vR*sin(pi*k)/(d*fD);

% TOTAL FORCE in recoil momentum/drive period
Ftot = @(x) (F-m*w^2*d*x)/(hbar*kL*fD);

% semiclassical transport: xdot = group velocity, hbar*kdot = force
dfdt = @(t,f) [vg(Vt(t),f(2));Ftot(f(1))];

%% Prepare Initial State
% phase space distribution
psd = @(x,k) exp(-x.^2/(2*(xwidth/d)^2))*exp(-k.^2/(2*kwidth^2));

% initial meshes (x labels rows, k labels columns)
xinit = linspace(-xend,xend,xpoints)'/d;       % initial x states
kinit = linspace(-kend,kend,kpoints);          % initial k states
pinit = psd(xinit,kinit);                      % probability grid
pinit = pinit/sum(pinit,'all');
dx = 2*xend/d/(xpoints-1);
dk = 2*kend/(kpoints-1);

% cell array to hold trajectories
xt = cell(xpoints,kpoints); kt = cell(xpoints,kpoints);

%% Integrate
tic
parfor idx = 1:xpoints*kpoints
    [ii,jj] = ind2sub([xpoints kpoints],idx);
    thisx = dx*(ii-(xpoints+1)/2); thisk = dk*(jj-(kpoints+1)/2);
    [~,temp] = ode45(@(t,f) dfdt(t,f),fD*[-tpre 0],[thisx thisk],opts);
    [~,Y] = ode45(@(t,f) dfdt(t,f),fD*time,[temp(end,1) temp(end,2)],opts);
    xt{idx} = Y(:,1); kt{idx} = Y(:,2); 
end
toc
xt = cell2mat(cellfun(@(x) reshape(x,1,1,[]),xt,'un',0));
kt = cell2mat(cellfun(@(x) reshape(x,1,1,[]),kt,'un',0));

%% Compute observables
xmean = squeeze(sum(repmat(pinit,1,1,length(time)).*xt,[1 2]));
xsigma = squeeze(sum(repmat(pinit,1,1,length(time)).*(xt.^2),[1,2]));
xsigma = sqrt(xsigma-xmean.^2);
xmean = xmean-xmean(1);

%% Plot
figure(57); clf;
subplot(211);
hold on;
plot(1E3*time,xmean,'m','linewidth',3);
xlabel('Time (ms)'); ylabel('Mean Position ($d$)','interpreter','latex');
subplot(212); hold on;
plot(1E3*time,xsigma,'m','linewidth',3);
xlabel('Time (ms)'); ylabel('Spatial Width ($d$)','interpreter','latex');

%% Drive Map
cycle = TD/tstep;
figure(357); clf;
% yline((F-h*fD/d)/(m*w^2*d),'k','linewidth',2); hold on;
[Xinit,Kinit] = meshgrid(xinit,kinit);
scatter(Kinit(:),Xinit(:),2,Kinit(:),'filled');
ylim([-500 500]); xlim([-1 1]);
ylabel('Position ($d$)','interpreter','latex'); xlabel('Quasiomentum ($k_L$)','interpreter','latex');
title('Curvature $f_D=53.56$ Hz, $T_B = 16.75$ ms, $\varphi=\pi/2$','interpreter','latex');
colormap(cmocean('balance'))
c = colorbar();
t=text(.45,.95,'');
for ii=1:tend/TD
%     clf
    hold on;
    ylim([-500 500]); xlim([-1 1]);
    xs = xt(:,:,ii)'; ks = mod(kt(:,:,ii)'+1,2)-1;
    scatter(ks(:),xs(:),2,Kinit(:),'filled');
    delete(t);
    t=text(.42,.95,['Period ' num2str(ii*cycle)],'fontname','times','fontsize',14,'units','normalized');
    drawnow;
    del = .05;
%     filename = 'curve_sparse.gif';
%     im = frame2im(getframe(gcf));
%     [imind,cm] = rgb2ind(im,256);
%     if ii == 1
%         imwrite(imind,cm,filename,'gif','Loopcount',Inf,'DelayTime',del);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
%     end
end