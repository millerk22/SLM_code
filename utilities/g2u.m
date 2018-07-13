function u = g2u(g,nhat)

	%[~,~,g] = unique(g); % Eliminate extraneous classes.

	if nargin == 1
		nhat = max(g);
	end

	N = length(g);
	u = zeros(N,nhat);
	for a=1:nhat;
		u((g==a),a) = 1;
	end
end
