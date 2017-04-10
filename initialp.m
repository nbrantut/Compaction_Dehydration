function [p,q,p1,q1] = initialp(pertp,  pm)
%INITIALP computes the initial total mean and shear stress across the
%layer, after a perturbation in pore pressure has been introduced. This
%ensures that the whole system remains on the yield curve.
%
%input:
%   pertp:  an array giving the pore pressure perturbation
%   pm:     the structure containing all the parameters.
%
%output:
%   p:   the array containing the total mean stress
%   q:   the array containing the shear stress
%   p1:  mean stress perturbation
%   q1:  shear stress perturbation

% compute normalised stress state:
pe0 = pm.pe0/abs(pm.sn0);
q0 = pm.q0/abs(pm.sn0);

%some useful notations
ps = pm.pstar0./pm.z0;
A = 3/4+pm.C^2;
B = q0*sqrt(3) - pm.C^2*(pm.b - 2*pe0 - ps - 2*pertp);
C = - pm.C^2*pertp.*(pm.b - 2*pe0 -ps - pertp);

D = B.^2 - 4*A.*C;

%get the perturbations:
p1 = (-B-sqrt(D))./(2*A);
q1 = sqrt(3)/2*p1;

%get the total values:
p = p1 + pm.p0./abs(pm.sn0);
q = q1 + q0;