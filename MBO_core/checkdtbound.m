%This checks the theoretical bounds from the MBO paper by Yves van Gennip to avoid trivial dynamics. These bounds are very far from sharp! I did modify the proofs to include the forcing terms, so the bounds I use here are slightly different. See my paper for details.
% dt is the time between thresholding, W is the adjacency matrix, and nhat is the number of classes to be constructed.

function dt = checkdtbound(dt,W,nhat,gamma,D,opts);

	%clean
	if min(size(D,1),size(D,2)) > 1
		D=diag(D)
	end

	N = size(W,1);
	
	% Weird case
	if N==2
		return;
	end

	k=sum(W,2);
	twom=sum(k);

	% Check that pinning will not automatically occur
	L=spdiags(k,0,N,N)-W;
	if N >  10
		M = @(u) L*u + 2*gamma/twom*k*(k'*u);
		opts.tol = 0.1;
		rho = eigs(M,N,1,'lm',opts);
	else
		M = L + 2*gamma/twom*k*k';
		[VM,DM] = eig(full(M));
		DM = diag(DM);
		rho = max(DM(DM>0));
	end
	t_rho=rho^(-1)*log(1 + .5*( nhat/N )^(0.5));

	% Check that steady state will not be achieved before thresholding. We assume that the volume of a cluster is not likely to substantially exceed N/nhat.
	lam_1 = D(1);
	lam_2 = D(2);
	tol = 10^(-12);
	if min(min(D)) < tol % D was not really supplied
		[~,D(1:2)] = my_eigs(W,gamma,k,twom,2,opts.thread,'M');
	end
	if lam_1 > tol
		t_t = 1/(2*lam_1) * log( 4*N / nhat );
	elseif lam_2 > 0 %the estimate from yves' paper works here
		if nhat > 2
			R_S = 1/nhat;
			t_t = 1/lam_2 * log( sqrt(N/nhat * (N-N/nhat)) / ( sqrt(N) * abs(R_S - 0.5) ) );
		else
			return; %Give up
		end
	else
		return;
	end

	%flag=1;
	%If a bound is violated, do something about it.
	flag=0;
	if t_rho > dt
		% disp('MBO thresholding too frequent. Pinning is guaranteed. t_rho')
		flag=1;
	end
	if dt > t_t
		% disp('Timestep too long for MBO. Trivial dynamics likely.')
		flag=1;
	end

	% Fix it if possible
	if flag==1
		if t_t <= t_rho
			% Not sure when this can occur
			disp('No way to prevent weird MBO behavior')
			% I do not expect this to occur
			disp('This is for unknown reasons')
			t_rho
			t_t
			lam_1
			lam_2
		else
			% disp('Adjusting timestep to avoid this issue')
			dt = sqrt(t_rho*t_t);
		end
	end
