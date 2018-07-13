Nfull=length(params.fulltag);
singletons=setdiff(1:Nfull,params.ind);
gfull(params.ind)=g;
gfull(singletons)=max(g)+(1:length(singletons));
g=gfull;
