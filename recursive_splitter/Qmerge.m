function [Cut,vol,Q] = Qmerge(Cut_old,vol_old,ahat,bhat,twom,gamma,opts)

	% Make sure bhat > ahat
	if bhat<ahat
		t = ahat;
		ahat = bhat;
		bhat = t;
	end

	nhat = size(Cut_old,2); %nhat is the size of the OLD partition

	ind = [1:bhat-1 bhat+1:nhat]; % Everything except bhat

	vol = vol_old(ind);
	vol(ahat) = vol(ahat) + vol_old(bhat);

	Cut = Cut_old(ind,ind);
	Cut(:,ahat) = Cut_old(ind,ahat)  + Cut_old(ind,bhat);
	Cut(ahat,:) = Cut_old(ind,ahat)' + Cut_old(ind,bhat)';
	Cut(ahat,ahat) = Cut_old(ahat,ahat) + Cut_old(bhat,bhat) + 2*Cut_old(ahat,bhat);

	% This seems bug-prone. Why not keep the canonical implementation in W_opt?
	tol = 10^(-7);
	vv = vol'*vol;
	iind = ((vol'*vol).*Cut) > tol;

	Q = sum(sum(Cut)) - sum(diag(Cut))  + (gamma/twom)*sum(vol.^2);

	nhat = size(Cut,2);

	if nargin < 4 || isfield(opts,'nhat') == false % Use AIC if no prior is provided.
		%warning('nhat was not supplied in Qmerge. Are you sure this is right?')
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
