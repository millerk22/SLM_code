function ic = Modularity_MBO(W,u0,gamma,V,D)
	[ic,dQ] = Modularity_MBO_sub(W,u0,gamma,V,D,sum(W,2),sum(sum(W)));
end
