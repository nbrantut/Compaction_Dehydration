function [hfig] = yieldsurfplot(sol,Im,pm,hfig)
%hfig = YIELDSURFPLOT(sol,Im,pm,hfig)
%
%YIELDSURFPLOT generates a plot of the trajectory, in the stress space, 
%of the point at index Im within the layer.
%
%input:
%   sol:  solution structure of the ODE
%   Im:   index of the point(s) to track
%   pm:   parameters structure
%   hfig: handle to figure object (optional)
%
%output:
%   hfig: handle to figure object.

% check input and make figure if needed
if nargin==3
    hfig = figure;
else
    figure(hfig);
end

% extract variables to ease manipulation
I = pm.I;
t   = sol.x;
xi  = sol.y(1:I,:);
pf  = sol.y(I+1:2*I,:);
eps = sol.y(2*I+1:3*I,:);
p   = sol.y(3*I+1:4*I,:);
z   = sol.y(4*I+1:5*I,:);
q   = sol.y(5*I+1:6*I,:);
sn  = sol.y(6*I+1:7*I,:);

%get the effective stress:
pe = p+pf;

Imid = (pm.I+1)/2;

%generate variables to plot initial and final yield curves
tabpe0 = linspace(-pm.pstar0/pm.z0, pm.b);
tabq0 = pm.C.*sqrt((pm.b - tabpe0).*(tabpe0 + pm.pstar0./pm.z0));
plot(tabpe0, tabq0, 'k--');
hold on;
tabpe = linspace(-pm.pstar0/z(Imid,end), pm.b);
tabq = pm.C.*sqrt((pm.b - tabpe).*(tabpe + pm.pstar0./z(Imid,end)));
plot(tabpe, tabq, 'k-');

%plot trajectory
plot(pe(Im,1), q(Im,1), 'ro',...
    'markersize',4,...
    'markeredgecolor','r');
plot(pe(Im,end), q(Im,end), 'ro',...
    'markersize',4,...
    'markeredgecolor','r');
plot(pe(Im,:)', q(Im,:)','r-');

axis equal;
set(gca, ...
    'xlim', 1.2*[tabpe0(1)  tabpe0(end)],...
    'ylim', [0  max(tabq0)*1.2]);
xlabel('({\itp}+{\itp}_f)/\sigma_n');
ylabel('\tau/\sigma_n');