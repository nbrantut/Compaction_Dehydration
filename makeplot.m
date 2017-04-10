if not(exist('sol','var') && exist('pm','var'))
    run ../code/mainq.m;
    close all;
end

figname = 'epsilon';
vr = eps;
varname = '\epsilon';%'{\itp}_f/\sigma_n';

ind = find(t<40,1,'last');

[~,ht1,ht2] = solutionplot(x,t,vr,varname,ind);
% set(ht1,'horizontalalignment','right');
% set(ht2,'horizontalalignment','right');


exportfig(figname,'xsize',13);