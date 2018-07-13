function [vol,Cut,Q] = Qmerge(vol_old,g_old,ig,ic_sub,k,Cut_old,A,opts)

	% Get volume
	nodes_ig = find(g_old==ig);
	g_new = ic2g(g_old,ic_sub,nodes_ig);

	groups_old = [1:ig-1 ig+1:max(g_old)];
	groups_new = [ig max(g_old)+1:max(g_new)];

	vol = zeros(1,max(g_new));
	vol(groups_old) = vol_old(groups_old);
	for i=groups_new
		vol(i) = sum(k(g_new == i));
	end

	% Get cut
	u_new = g2u(g_new);
	nodes_old = find(sum(u_new(:,groups_old),2));
	nodes_new = find(sum(u_new(:,groups_new),2));

	Cut = zeros(max(g_new));
	Cut(groups_old,groups_old) = Cut_old(groups_old,groups_old);

	us = u_new(nodes_new,groups_new);
	Cut(groups_new,groups_new) = us'*A(nodes_new,nodes_new)*us;
	
	up = u_new(nodes_old,groups_old);
	Cut(groups_new,groups_old) = us'*A(nodes_new,nodes_old)*up;

	Cut = max(Cut,Cut');

	nhat = size(Cut,2);
	u = g2u(g_new,nhat); % the nhat is important in this case

	k = full(sum(A,2));
	twom = sum(k);
	m = twom/2;

	Q = sum(sum(Cut)) - sum(diag(Cut)) + (gamma/twom)*sum(sum(vol'*vol));

	if isfield(opts,'nhat') == false % Use AIC if no prior is provided.
		Q = Q + nhat^2;
		return;
	else
		nhat_opt = opts.nhat;
	end

	if isfield(opts,'lam') == false
		lam = 0.1;
	else
		lam = opts.lam;
	end

	Q = Q + lam*abs(Q)*(nhat - nhat_opt)^2;
end

function g_temp = ic2g(g,ic_sub,ind)
	g_temp = g;
	for j=2:max(ic_sub)
		tempInd_temp=(ic_sub==j);              
		g_temp(ind(tempInd_temp))= max(g)+j-1;              
	end
end
