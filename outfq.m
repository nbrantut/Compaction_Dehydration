function status = outfq(t,y,flag,params)
%OUTFQ is the output fonction used while solving the ODE. It plots the
% current yield curve at the middle of the layer, and points to the current
% stress state there. In order to check consistency, the shear stress
% computed from the geometrical constrain is also plotted, and should
% overlap with the one obtained from the ODE solution.

if strcmp(flag,'init')
    tabpe0 = linspace(-params.pstar0/params.z0, params.b);
    tabq0 = params.C.*sqrt((params.b - tabpe0).*(tabpe0 + params.pstar0./params.z0));
    plot(tabpe0, tabq0, 'k--');
    axis equal;
    set(gca, ...
        'xlim', 1.2*[tabpe0(1)  tabpe0(end)],...
        'ylim', [0  max(tabq0)*1.2]);
    xlabel('({\itp}+{\itp}_f)/\sigma_n');
    ylabel('{\itq}/\sigma_n');
    hold on;
    
end
 
if not(strcmp(flag,'done'))
    
    %assign variable names for simplicity
    I = params.I;

    pf  = y(I+1:2*I);

    p   = y(3*I+1:4*I);
    z   = y(4*I+1:5*I);
    q   = y(5*I+1:6*I);

    %get the shear stress:
    q2 = -sqrt(3)/2*(params.sn - p);
    
    %and the effective stress:
    pe = p+pf;

    Imid = (params.I+1)/2;

    tabpe = linspace(-params.pstar0/z(Imid), params.b);
    tabq = params.C.*sqrt((params.b - tabpe).*(tabpe + params.pstar0./z(Imid)));
    h2=plot(tabpe, tabq, 'k-');
    plot(pe(Imid), q(Imid), 'ro',...
        'markerfacecolor','r',...
        'markersize',5)
    h3=plot(pe(Imid), q2(Imid), 'b+');
    title(['{\itt} = '  num2str(t)]);
    
    drawnow;
    
    delete([h2 h3]);
    
else
    close;
    
end

status = 0;