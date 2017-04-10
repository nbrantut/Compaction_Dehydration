function [hfig,ht1,ht2] = solutionplot(x,t,vr,varname,ind)

hfig = figure;

Imid = (length(x)+1)/2;

subplot 121

plot(x,vr(:,1));
hold all;

cm = colormap;
N = length(cm);
step = 25;
Nt = length(t);
tabk = 1:step:Nt;
Nk = length(tabk);

for k=1:Nk
    ic = ceil(N/Nk*k);
    plot(x,vr(:,tabk(k)), 'color',cm(ic,:));
end

bounds=[min(min(vr)) max(max(vr))];
rg = diff(bounds);
ext = rg/20;

xlim([x(1) x(end)])
ylim(bounds + [-ext ext]);
xlabel('{\ity}/{\itL}');
ylabel(varname);

set(gca, 'position', get(gca,'position').*[1 1.5 1 1])

subplot 122

plot(t,vr(Imid,:),'k');
hold on;
plot(t, vr(1,:),'k');
xlabel('{\itt}\times{\itr}_0');
ylim(bounds + [-ext ext]);

set(gca, 'yticklabel', [],...
    'position', get(gca,'position').*[0.9 1.5 1 1])

ht1=text(t(ind), vr(Imid,ind), '{\ity}=0',...
    'verticalalignment','bottom');
ht2=text(t(ind), vr(1,ind), '{\ity}=\pm {\itL}/2',...
    'verticalalignment','bottom');