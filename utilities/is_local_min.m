function [flag,i,l] = is_local_min(g,A,gamma)
	flag=true;

	N = length(g);
	[~,~,g] = unique(g);
	nhat = max(g);

	k = full(sum(A,2));
	twom = sum(k);

	u=g2u(g);
	Cut = u'*A*u;

	vol = k'*u;

	Qold = sum(sum(Cut)) - sum(diag(Cut)) + gamma/twom*sum(vol.^2);

	for i=1:N
		for l=1:nhat
			if g(i) ~= l

				gtemp = g;
				gtemp(i)=l;
				u=g2u(gtemp);
				Cut = u'*A*u;

				vol = k'*u;

				Q = sum(sum(Cut)) - sum(diag(Cut)) + gamma/twom*sum(vol.^2);

				tol = 10^(-9);
				if Qold > Q + tol
					flag=false;
					return;
				end
			end
		end
	end

end
