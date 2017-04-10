%% Numerical computation of dehydration/compaction instability
% 
% Here we compute numerically the evolution of pore pressure, reaction
% progress, stress and strain for an antigorite layer in uniaxial
% deformation.
%
% Specifically, this code solves the following set of coupled PDE/ODEs:
% 
% $$
% \frac{\partial\xi}{\partial t} = s(1-\xi)|p_f/p_{eq}-1|^{n_r},
% $$
% 
% $$
% \frac{\partial p_f}{\partial t} =
% \frac{Mk/\eta}{1+\beta^*M}\frac{\partial^2 p_f}{\partial y^2} + \frac{M
% m_d^0 (1/\rho_f +\Delta_rV_s) - f'X}{1+\beta^*M}\frac{\partial
% \xi}{\partial t},
% $$
% 
% $$
% \frac{\partial \epsilon}{\partial t} = \frac{1}{M}\left(\frac{\partial 
% p_f}{\partial t} + f'X\frac{\partial\xi}{\partial t}\right),
% $$
% 
% $$
% \frac{\partial p}{\partial t} =
% \frac{1}{\sqrt{3}/2+\beta-f'/K}\left((f'/K-\beta)\frac{\partial
% p_f}{\partial t} - f'\frac{\partial\epsilon}{\partial t} +
% f'm_d^0\Delta_rV_s\frac{\partial\xi}{\partial t}\right),
% $$
% 
% $$
% \frac{\partial\zeta}{\partial t} =
% -m_d^0\Delta_rV_s\frac{\partial\xi}{\partial t} +
% \frac{\partial\epsilon}{\partial t} - \frac{1}{K}\left( \frac{\partial
% p}{\partial t}+\frac{\partial p_f}{\partial t}\right),
% $$
% 
% $$
% \frac{\partial \tau}{\partial t} = -\beta\left(\frac{\partial p}{\partial
% t}+ \frac{\partial p_f}{\partial t}\right) -
% f'\frac{\partial\xi}{\partial t}.
% $$
% 
% In addition, a consistency check is performed by computing the total
% normal stress:
% 
% $$
% \frac{\partial\sigma_n}{\partial t} = \frac{\partial p}{\partial t} -
% \frac{2}{\sqrt{3}}\frac{\partial\tau}{\partial t},
% $$
% 
% and verifying a posteriori that it remains constant throughout space and
% time.
%
% All the equations and parameters are normalised, using $|\sigma_n|$ as 
% the stress and pressure scale,  $1/r_0$ as the time scale, and the layer 
% width as the space scale.

%% Set parameters
%
% Here we use the function |parameters| to get all the parameter values and
% functions required for the computation.

% initial conditions:
z0 = 0.03;     %porosity
peq = 3e9;     %equilibrium pressure
p0 = peq/0.828; %mean stress

% fetch parameters:
pm = parameters(...
    'z0',  z0,...
    'peq', peq,...
    'p0',  p0,...
    'L', 100);

%% Check initial conditions and stress state
%
% Here we basically display some information about the initial (normalised)
% stress state and parameter values.

% display some information about the initial conditions:
disp('Initial stress conditions:');
disp(['  q  = '  num2str(pm.q0/1e9)   ' GPa']);
disp(['  p  = '  num2str(pm.p0/1e9)   ' GPa']);
disp(['  sn = '  num2str(pm.sn0/1e9)  ' GPa']);
disp(['Dilatancy factor: \beta = '  num2str(pm.B(-pm.pe0/pm.sn0,-pm.q0/pm.sn0,pm.z0))]);

% plot the initial yield surface and where we are:
tabpe = linspace(-pm.pstar0/pm.z0, pm.b);
tabq = pm.C.*sqrt((pm.b - tabpe).*(tabpe + pm.pstar0./pm.z0));
plot(tabpe, tabq, 'k-');
hold on;
plot(-pm.pe0./pm.sn0, -pm.q0./pm.sn0, 'ro');
axis equal;
set(gca, ...
    'xlim', 1.2*[tabpe(1)  tabpe(end)],...
    'ylim', [0  max(tabq)*1.2]);
xlabel('({\itp}+{\itp}_f)/\sigma_n');
ylabel('{\itq}/\sigma_n');

%% Prepare computation and solver options
%
% This is to prepare all the required variables to input in the ODE solver.

%get I to ease notations:
I = pm.I;
x = pm.x;
pertp = 1e-6*(1+cos(2*pi*x'));

%initialise solution
xi0  = zeros(I,1);
pf0  = pm.peq + pertp;
eps0 = zeros(I,1);
z0   = pm.z0 + zeros(I,1);
sn0  = pm.sn + zeros(I,1);
[p0,q0,~,~] = initialp(pertp,  pm);

y0 = [xi0; pf0; eps0; p0; z0; q0; sn0];

%time interval:
tspan = [0 10000];

%ode options:
options = odeset(...
    'AbsTol', 1e-8,...
    'RelTol', 1e-6,...
    'OutputFcn', @(t,y,flag) outfq(t,y,flag,pm));
    

%% Solve and extract results
%
% The solution is obtained using the stiff ode solver |ode15s|. Tests with
% |ode23s| or |ode23t| turned out to be very slow.

sol = ode15s(@(t,y) fqp(y, pm), tspan, y0, options);

% extract variable to ease manipulation:
t   = sol.x;
xi  = sol.y(1:I,:);
pf  = sol.y(I+1:2*I,:);
eps = sol.y(2*I+1:3*I,:);
p   = sol.y(3*I+1:4*I,:);
z   = sol.y(4*I+1:5*I,:);
q   = sol.y(5*I+1:6*I,:);
sn  = sol.y(6*I+1:7*I,:);

% get the index of the central node:
Imid = (pm.I+1)/2;

%% Plot results
%
% Here we compute selected plots showing the key results.

%trajectory in the stress space 
yieldsurfplot(sol,(1:I),pm);



