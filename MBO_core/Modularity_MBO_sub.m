% This is the version to be used in the recursive scheme. It differs mostly in small details.
% It is important to note that k is the degree of the node in the original (large) graph, rather than the subgraph we are working in.
% I am now modifying the code so that it accepts and works with the spectral data for M = L + gamma/m kk^T. This is true only for the balanced tv version.
function [groups,dQ] = Modularity_MBO_sub(W_sub, u0, gamma,V,D,k,twom)

	%housekeeping
	if(min(size(D))>1)
		D=diag(D);
	end

	% parameters
	iter =1000;
	twom_sub = sum(k);
	N = length(k);
	nhat=size(u0,2);

	% initialization
	u = u0;
	[~,PreIndex] = max(u,[],2);

	dt_thresh = 12*.02;% How frequently to threshold
	dt_split = .02;% Timestep for diffusion and forcing
	s = dt_thresh/dt_split; % How many timesteps to take before thresholding.
	dt = dt_split;
	%s=5;
	%dt=1;

	Denom = 1 + dt*D;
	a = V'*u;

	% main iteration
	for j = 1:iter
		for l=1:s
			meanu = (k'*u)/twom_sub;
			a = V'*u;
			temp = 2*gamma*twom_sub*(spdiags(k,0,N,N)*(u-ones(N,1)*meanu))/twom;
			b = V'*temp;
			a = diag((Denom).^(-1))*(a +dt*b);
			u = V*a;                  
		end

		%MBO thresholding:
		[~,Index] = max(u,[],2);
		u = sparse(1:N,Index,ones(N,1),N,nhat);

		%Stopping Criterion -- try more frequent thresholding if you get stuck
		if norm(PreIndex-Index)/norm(Index)<.001 %|| mod(j,50)==0
			%if dt< 10^(-3);%10^(-2) works well here
				break;
			%else
				%dt =.4*dt;%dt=.4*dt works well here. .2 also works.
			%end
		else
			PreIndex = Index;
		end;

	end

	if j==iter
		disp('MBO did not converge'); 
	end

	[~,~,groups]=unique(Index);

	Q_div= sum(sum(u.*(W_sub*u-(gamma/twom)*k*(k'*u))))/twom;
	Q_nodiv = (sum(sum(W_sub))-gamma*sum(k)^2/twom)/twom;
	dQ=Q_div-Q_nodiv;
