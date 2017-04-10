function param = parameters(varargin)
%PARAMETERS generates a structure containing all the required input
%parameters to solve the compaction/dehydration problem.
%
%It accepts input in the form of 'param', value pairs.

disp('Assigning parameters...');

p = inputParser;

p.addParamValue('I',64+1, @isnumeric);
p.addParamValue('L',20, @isnumeric);
p.addParamValue('atg',1, @isnumeric);
p.addParamValue('K',67.9e9, @isnumeric);
p.addParamValue('G',38.5e9, @isnumeric);
p.addParamValue('r0',1.47e-6, @isnumeric);
p.addParamValue('nr',1, @isnumeric);
p.addParamValue('k0',1e-20, @isnumeric);
p.addParamValue('C',0.8, @isnumeric);
p.addParamValue('b', 5e6, @isnumeric);
p.addParamValue('pstar0', 19e6, @isnumeric);
p.addParamValue('p0',3e9/0.9, @isnumeric);
p.addParamValue('h',0,@isnumeric);
p.addParamValue('z0',0.04, @isnumeric);
p.addParamValue('peq',3e9, @isnumeric);

%parse input parameters:
p.parse(varargin{:});
param = p.Results;

%get stress conditions:
disp('Get stress conditions...');
param.p0 = -param.p0;   %minus sign here so that compression <0
param.pe0 = param.p0 + param.peq;

%check if we are good:
if abs(param.pe0)>=abs(param.pstar0./param.z0)
    disp('Effective stress is too large ! Change p or z0.');
    disp('Abort !');
    return
end

param.q0 = param.C.*sqrt((param.b - param.pe0)...
    .*(param.pe0 + param.pstar0./param.z0));
param.sn0 = param.p0 - 2/sqrt(3)*param.q0;

%get properties that depend on pressure:
disp('Run thermodynamic calculations...');
run ../figures/volumes;
disp('Extract PT-dependent properties...');
bf = interp1(r2.pb(:,2), r2.bf, param.peq/1e9);
eta = interp1(r2.pb(:,2), r2.etapb, param.peq/1e9);
rhof = interp1(r2.pb(:,2), 1./r2.Vh2o, param.peq/1e9);
DrVs = interp1(r2.pb(:,2), r2.DrVspb, param.peq/1e9);
md0tot = interp1(r2.pb(:,2), r2.md0pb, param.peq/1e9);
param.betaf = bf;
param.eta = eta;
param.rhof = rhof;
param.DrVs = DrVs;
param.md0 = param.atg.*md0tot;

%normalise all stresses by |\sigma_n|, times by 1/r0 and lengths by L:
disp('Normalise parameters...');
param.K = param.K/abs(param.sn0);
param.G = param.G/abs(param.sn0);
param.peq = param.peq/abs(param.sn0);
param.betaf = param.betaf*abs(param.sn0);
param.b = param.b/abs(param.sn0);
param.pstar0 = param.pstar0/abs(param.sn0);
param.h = param.h/abs(param.sn0);
param.sn = param.sn0/abs(param.sn0);

param.eta = param.eta/abs(param.sn0)*param.r0;
param.k0 = param.k0/param.L^2;

param.x = linspace(-1/2,1/2,param.I);
param.Dx = param.x(2)-param.x(1);


disp('Generate functions...');
param.k = @(phi) param.k0*((1-param.z0)./(1-phi)).^2.*(phi./param.z0).^3;
param.betas = @(phi) (1-phi)./param.K + phi.*param.betaf;

param.f = @(pe,q,z) q - param.C.*sqrt(...
    (param.b - pe).*(pe + param.pstar0./z)...
    );

param.B = @(pe,q,z) -(param.C.*(param.b - 2.*pe - param.pstar0./z))./...
    ( 2 * sqrt((param.b - pe).*(pe + param.pstar0./z)) );

param.dfdz = @(pe,q,z) (param.C.*(param.b - pe).*param.pstar0)./...
    (2*sqrt((param.b - pe).*(pe + param.pstar0./z)).*z.^2);

disp('Done!');
