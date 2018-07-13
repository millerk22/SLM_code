function Q = Qcut(g,A,gamma,Cut,vol)

	nhat = size(gamma,2);
	u = g2u(g,nhat); % the nhat is important in this case

	k = full(sum(A,2));
	twom = sum(k);

	if nargin == 3
		Cut = u'*A*u;
		vol = k'*u;
	end

	Q = sum(sum(Cut)) - sum(diag(Cut)) + (gamma/twom)*sum(vol.^2);

end
