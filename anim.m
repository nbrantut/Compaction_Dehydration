function hf = anim(x,t,vr,step,time,varname)

if nargin==3
    step=1;
    time=1e-3;
    varname='';
elseif nargin==4
    time=1e-3;
    varname='';
elseif nargin==5
    varname='';
end

hf = figure;
plot(x,vr(:,1));
xlim([x(1) x(end)])
ylim([min(min(vr)) max(max(vr))]);
xlabel('{\ity}/{\itL}');
ylabel(varname);
set(gca,'nextplot','replacechildren');
for k=1:step:length(t)
    plot(x,vr(:,k));
    title(['{\itt}/\tau = '  num2str(t(k))]);
    drawnow;
    pause(time);
end