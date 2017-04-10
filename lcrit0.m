function out = lcrit0(params)
%LCRIT0 computes the critical (normalised) wavelength for strain
%localisation.
% usage: out = LCRIT0(params)
%
%input:
%   params:  parameters structure
%
%output:
%   out:     critical wavelength

% assign variable names for simplicity
eta = params.eta;
md0 = params.md0;
rhof = params.rhof;
DrVs = params.DrVs;
peq = params.peq;
z = params.z0;
q = params.q0./abs(params.sn0);

% and the effective stress:
pe = params.pe0./abs(params.sn0);

% assign variable parameters
bs = params.betas(z);
koe = params.k(z)/eta;
B = params.B(pe,q,z);
mu = B;
h = params.h;
G = params.G;
K = params.K;
fp = params.dfdz(pe,q,z);

% compute some useful lumped parameters:
M = G*K./(h-fp.*B+G+B.*mu.*K) .*((1+2/sqrt(3)*B).*(1+2/sqrt(3)*mu) ...
    + (h-fp.*B).*(1/G + 4/3/K));
X = -md0*DrVs*(B*K - 2*G/sqrt(3))./(h-fp.*B + G +B.*mu.*K);

out = 2*pi*sqrt(-pe*koe./(md0*(1/rhof+DrVs)-fp.*X./M));