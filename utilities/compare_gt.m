load pnnl_100_louvain_gamma_10
sig=csvread('100K_ind.csv');
[~,~,g]=unique(g(sig));
g=hist(g,max(g));
g=sort(g,'descend');
