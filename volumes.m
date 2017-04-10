%% Computation of mineral and fluid properties along the phase boundary
% 
% This script is used to determine all the mineral and fluid volume,
% compressbility and the fluid viscosity along the phase boundary between
% antigorite (atg) and
% # talc (tlc) and forsterite (fo),
% # enstatite (ens) and forsterite,
% # phase A (phA) and enstatite.
%
% It outputs three structures |r1|, |r2| and |r3| containing the results
% for each reaction.
%
% All the data, except when explicitly mentioned otherwise, are from the 
% table in Hacker et al., JGR 2003.
%
% The water properties are computed using the IAPWS formulation,
% implemented in the Matlab Package water95.

%% Mineral parameters

%antigorite
atg.V0 = 1754.7e-6;
atg.M = 4536e-3;
atg.a = 4.7e-5;
atg.K = 6.79e10; %this value is not the same as in Hacker, this is 
                 %for consistency with the Bezacier values.
atg.dKdP = 2.77;

%forsterite
fo.V0 = 43.7e-6;
fo.M = 140.7e-3;
fo.a = 6.1e-5;
fo.K = 12.7e10;
fo.dKdP = 5.37;

%enstatite
ens.V0 = 62.6e-6;
ens.M = 200.8e-3;
ens.a = 5.1e-5;
ens.K = 10.6e10;
ens.dKdP =  8.5;

%talc
tlc.V0 = 136.4e-6;
tlc.M = 379.7e-3;
tlc.a = 3.7e-5;
tlc.K = 4.16e10;
tlc.dKdP = 6.5;

%phase A
phA.V0 = 154.4e-6;
phA.M = 456.3e-3;
phA.a = 8.3e-5;
phA.K = 9.74e10;
phA.dKdP = 6.0;

%water
h2o.M = 18e-3;

%% Functions to compute density

T0 = 25; %Celcius

rhoTnorm = @(a,T) exp(-a*(T-T0 -20*(sqrt(T)-sqrt(T0)) ));

strain = @(f,dKdP) 3*f.*(1+2*f).^(5/2) .*(1-2*(3-.75*dKdP)*f + ...
    (f.^2 /6).*(4*(3-.75*dKdP)*(4-3*dKdP) + 5*(3*dKdP-5)));

ff = @(P,K,dKdP) fzero(@(f) P/K - strain(f,dKdP) ,0);

rhoPnorm = @(P,K,dKdP) (1 + 2*ff(P,K,dKdP)).^(2/3);

Vm = @(P,T,Vm0,a,K,dKdP) Vm0./(rhoTnorm(a,T) .* rhoPnorm(P,K,dKdP));

betaf = @(p,T) (1/density(p,T))*(density(p+1,T)-density(p-1,T))/2;

%% Reactions

% reaction 1: to tlc, fo, water
r1.md0  = @(Vm_atg) 27.*h2o.M./Vm_atg;
r1.DrVs = @(Vm_tlc,Vm_fo,Vm_atg) (4*Vm_tlc + 18*Vm_fo - Vm_atg)./(27*h2o.M);

% reaction 2: to fo, ens, water
r2.md0  = @(Vm_atg) 31.*h2o.M./Vm_atg;
r2.DrVs = @(Vm_ens,Vm_fo,Vm_atg) (10*Vm_ens + 14*Vm_fo - Vm_atg)./(31*h2o.M);

% reaction 3: to phA, ens, water
r3.md0  = @(Vm_atg) (113/5).*h2o.M./Vm_atg;
r3.DrVs = @(Vm_ens,Vm_phA,Vm_atg) ((71/5)*Vm_ens + (14/5)*Vm_phA - Vm_atg)./((113/5)*h2o.M);

%% Phase boundaries
% Here we use the script |phaseboundaries.m| which contains the table of 
% pressure and temperature coordinates of the phase boundary.

phaseboundaries;

pp = linspace(dat(1,2),8,500);

r1.pb(:,2) = pp(pp<1.61);
r1.pb(:,1) = interp1(dat(1:16,2),dat(1:16,1),r1.pb(:,2));

r2.pb(:,2) = pp(pp>=1.61 & pp<5.92);
r2.pb(:,1) = interp1(dat(17:38,2),dat(17:38,1),r2.pb(:,2),...
    'linear','extrap');

r3.pb(:,2) = pp(pp>=5.92);
r3.pb(:,1) = interp1(dat(39:end,2),dat(39:end,1),r3.pb(:,2),...
    'linear','extrap');

for i=1:length(r1.pb)
    r1.Vh2o(i) = 1./density(r1.pb(i,2)*1e9,r1.pb(i,1)+273.15);
    r1.DrVspb(i) = r1.DrVs(Vm(r1.pb(i,2)*1e9,r1.pb(i,1),...
        tlc.V0,tlc.a,tlc.K,tlc.dKdP),...
        Vm(r1.pb(i,2)*1e9,r1.pb(i,1),fo.V0,fo.a,fo.K,fo.dKdP),...
        Vm(r1.pb(i,2)*1e9,r1.pb(i,1),atg.V0,atg.a,atg.K,atg.dKdP));
    r1.md0pb(i) = r1.md0(Vm(r1.pb(i,2)*1e9,r1.pb(i,1),...
        atg.V0,atg.a,atg.K,atg.dKdP));
    r1.etapb(i) = visc_w(r1.pb(i,1)+273.15,1/r1.Vh2o(i));
    r1.bf(i) = betaf(r1.pb(i,2)*1e9, r1.pb(i,1)+273.15);
end

for i=1:length(r2.pb)
    r2.Vh2o(i) = 1./density(r2.pb(i,2)*1e9,r2.pb(i,1)+273.15);
    r2.DrVspb(i) = r2.DrVs(Vm(r2.pb(i,2)*1e9,r2.pb(i,1),...
        ens.V0,ens.a,ens.K,ens.dKdP),...
        Vm(r2.pb(i,2)*1e9,r2.pb(i,1),fo.V0,fo.a,fo.K,fo.dKdP),...
        Vm(r2.pb(i,2)*1e9,r2.pb(i,1),atg.V0,atg.a,atg.K,atg.dKdP));
    r2.md0pb(i) = r2.md0(Vm(r2.pb(i,2)*1e9,r2.pb(i,1),...
        atg.V0,atg.a,atg.K,atg.dKdP));
    r2.etapb(i) = visc_w(r2.pb(i,1)+273.15,1/r2.Vh2o(i));
    r2.bf(i) = betaf(r2.pb(i,2)*1e9, r2.pb(i,1)+273.15);
end

for i=1:length(r3.pb)
    r3.Vh2o(i) = 1./density(r3.pb(i,2)*1e9,r3.pb(i,1)+273.15);
    r3.DrVspb(i) = r3.DrVs(Vm(r3.pb(i,2)*1e9,r3.pb(i,1),...
        ens.V0,ens.a,ens.K,ens.dKdP),...
        Vm(r3.pb(i,2)*1e9,r3.pb(i,1),phA.V0,phA.a,phA.K,phA.dKdP),...
        Vm(r3.pb(i,2)*1e9,r3.pb(i,1),atg.V0,atg.a,atg.K,atg.dKdP));
    r3.md0pb(i) = r3.md0(Vm(r3.pb(i,2)*1e9,r3.pb(i,1),...
        atg.V0,atg.a,atg.K,atg.dKdP));
    r3.etapb(i) = visc_w(r3.pb(i,1)+273.15,1/r3.Vh2o(i));
    r3.bf(i) = betaf(r3.pb(i,2)*1e9, r3.pb(i,1)+273.15);
end
