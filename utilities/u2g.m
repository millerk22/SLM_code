function g = u2g(u)

	N = size(u,1);
	g = zeros(N,1);
	nhat = size(u,2);

	for a=1:nhat
		g(u(:,a)==1) = a;
	end

	%[~,~,g] = unique(g); % Eliminate extraneous classes.

end
