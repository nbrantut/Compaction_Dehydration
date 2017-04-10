function dy = fq(y, params)
%FQ is the function f(t,y) such that y' = f(t,y) is the system of ODE
%we want to solve.
%
%input:
%   y:   column vector containing, in that order, the reaction progress,
%   pore pressure, strain, mean stress, porosity, shear stress, and normal
%   stress.
%   params: structure (output from the function 'parameters') containing
%   all the nondimensional parameters and functions.
%
%ouput:
%   dy:  column vector containing the derivatives y'.

% Extract current variables and parameters:

% assign variable names for simplicity
I = params.I;
Dx = params.Dx;
eta = params.eta;
md0 = params.md0;
rhof = params.rhof;
DrVs = params.DrVs;
peq = params.peq;
nr = params.nr;

xi  = y(1:I);
pf  = y(I+1:2*I);
p   = y(3*I+1:4*I);
z   = y(4*I+1:5*I);
q   = y(5*I+1:6*I);

% and the effective stress:
pe = p+pf;

% assign variable parameters
bs = params.betas(z);
koe = params.k(z)/eta;
B = params.B(pe,q,z);
mu = B;
h = params.h;
G = params.G;
K = params.K;
fp = params.dfdz(pe,q,z);

% Safety check:
if any(imag(B))
    disp('B is complex ! Wrong!')
    return
end

% compute some useful lumped parameters:
M = G*K./(h-fp.*B+G+B.*mu.*K) .*((1+2/sqrt(3)*B).*(1+2/sqrt(3)*mu) ...
    + (h-fp.*B).*(1/G + 4/3/K));
X = -md0*DrVs*(B*K - 2*G/sqrt(3))./(h-fp.*B + G +B.*mu.*K);
denom = 1 + bs.*M;

% Derivatives

% reaction rate
dxi = sign(pf./peq-1).*(1-xi).*abs(pf./peq-1).^nr;

% pore pressure rate (using periodic BC)
dpf(1,1) = (M(1)./denom(1)).*(1/Dx).*(...
    0.5*(koe(2)+koe(1)).*(pf(2)-pf(1))/Dx ...
    -0.5*(koe(1)+koe(1)).*(pf(1)-peq)/Dx) ...
    + ( M(1)*md0*(1/rhof+DrVs)  - fp(1).*X(1))/denom(1).*dxi(1);

dpf(I,1) = (M(I)./denom(I)).*(1/Dx).*(...
    0.5*(koe(I)+koe(I)).*(peq-pf(I))/Dx ...
    -0.5*(koe(I)+koe(I-1)).*(pf(I)-pf(I-1))/Dx) ...
    + ( M(I)*md0*(1/rhof+DrVs)  - fp(I).*X(I))/denom(I).*dxi(I);

dpf(2:I-1,1) = (M(2:I-1)./denom(2:I-1)).*(1/Dx).*(...
    0.5*(koe(3:I)+koe(2:I-1)).*(pf(3:I)-pf(2:I-1))/Dx ...
    -0.5*(koe(2:I-1)+koe(1:I-2)).*(pf(2:I-1)-pf(1:I-2))/Dx) ...
    + ( M(2:I-1)*md0*(1/rhof+DrVs)  - fp(2:I-1).*X(2:I-1))./denom(2:I-1).*dxi(2:I-1);

% strain rate
deps = 1./M .*(dpf + fp.*X.*dxi);

% total mean stress rate
dp = 1./(sqrt(3)/2+B-fp./K).*((fp./K-B).*dpf -fp.*deps +fp*md0*DrVs.*dxi);

% inelastic porosity change
dz = -md0*DrVs*dxi + deps  - (dp+dpf)./K;

% shear stress rate
dq = -B.*(dp+dpf) - fp.*dz;

% normal stress
dsn = dp-2/sqrt(3)*dq;

% output vector of derivatives:
dy = [dxi;
    dpf;
    deps;
    dp;
    dz
    dq
    dsn];

end