% Based off of code writen by Hui Hu for the modularity MBO paper. Modified (a lot) by Zach Boyd.
% 
% A:  weight matrix
% gamma: resolution parameter
% n_given: max number of clusters at each sub-partition step. (One can choose n=10)
% Output g: group assignments


function [g,info] = recursive_splitter(A,gamma,n_given, opts)
	global verbose

	N = size(A,2);
	g = ones(N,1);

	info.nhat_max = 1;

	tol = 0.01;
	sims = -1;
	i = 1;
	cmin = 1;
	to_cont = true;
	to_split = true;
	while to_cont == true

		g_old = g;

		to_split = ones(max(g),1)' == 1;
		g = split(A,gamma,n_given,g,opts,to_split);

		g_old_ = g;

		changes(i) = nmiMeas(g,g_old);
		i = i+1;

		sims = nnz(g-g_old_);
		if verbose; disp(sprintf('sims %f',sims)); end

		if sims == 0
			to_cont = false;
		end

	end

	if length(changes) > 1
		changes(2:end);
	end

end
