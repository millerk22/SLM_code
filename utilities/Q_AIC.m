function Q = Q_AIC(g,A,W,opts,Cut,vol)

	% Note that if we are working with a subgraph, all 5 arguments must be specified, or the degree will be wrong!

	if nargin < 5
		Q = Qcut(g,A,W);
	else
		Q = Qcut(g,A,W,Cut,vol);
	end

end
