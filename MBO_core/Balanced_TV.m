% k is the degree of the node in the original (large) graph, rather than the subgraph we are working in.
function [groups,dQ] = Balanced_TV(W_sub, u0, gamma,V,D,k,twom,opts)
	pseudospectral = opts.pseudospectral;
	ssl = opts.ssl;

	%housekeeping
	if(min(size(D))>1)
		D=diag(D);
	end
	k=full(k);
	twom=full(twom);
	u0=full(u0);
	isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

	% parameters
	iter =1000;
	twom_sub = sum(k);
	N = length(k);
	nhat=size(u0,2);

	% initialization
	u = u0;
	[~,PreIndex] = max(u,[],2);

	for i=1:nhat % make sure each group is represented
		u(i,:)=0;
		u(i,i)=1;
	end

	% Makes sure the MBO timestep is reasonable.
	dt = .40; % default
	dt = checkdtbound(dt,W_sub,nhat,gamma,D,opts);

	if pseudospectral
		%disp('selecting pseudospectral solver')
	else
		%disp('using cg')
		L = spdiags(k,0,N,N) - W_sub;
	end
	if ssl.use_me
		%disp('using ssl');
	else
		%disp('omitting ssl');
	end

	% main iteration
	for j = 1:iter

		if ssl.use_me
			dti=min(1/(2*ssl.lam),dt);
			niter = ceil(dt/dti)
		else
			dti=dt;
			niter=1;
		end

		if pseudospectral
			DenomCV = diag(exp(-dti*D));
		else
			MM = @(x) reshape(reshape(x,[N,nhat]) + dti*L*reshape(x,[N,nhat]) + dti*gamma/twom*k*(k'*reshape(x,[N,nhat])),[N*nhat,1]);
		end

		for ii=1:niter

			if pseudospectral
				a = V'*u;
				a = DenomCV*a;
				u = V*a;
			else
				for zzz=1:10
					%disp('inner-loop timestep not theoretically justified');
					MM = @(x) reshape(reshape(x,[N,nhat]) + 0.1*dti*L*reshape(x,[N,nhat]) + 0.1*dti*gamma/twom*k*(k'*reshape(x,[N,nhat])),[N*nhat,1]);
					u=reshape(u,[N*nhat,1]);
					ttol=10^(-5);
					mmaxit=100;
					if isOctave
					  [u,flag,~,~,~,~] = pcg(MM,u,ttol,mmaxit);
					else
					  [u,flag,~,~,~] = pcg(MM,u,ttol,mmaxit);
					end
					if flag
						disp('pcg not converge')
					end
					u=reshape(u,[N,nhat]);
				end
			end

			if ssl.use_me
				u(ssl.ind,:) = u(ssl.ind,:) - dti*ssl.lam*(u(ssl.ind,:) - ssl.u(ssl.ind,:));
			end
		end

		%MBO thresholding:
		[~,Index] = max(u,[],2);
		u = sparse(1:N,Index,1,N,nhat);

		%Stopping Criterion -- try more frequent thresholding if you get stuck
		if norm(PreIndex-Index)/norm(Index)<.001 || mod(j,100)==0
			if dt< 10^(-3);%10^(-2) works well here
				break;
			else
				dt =.4*dt;%dt=.4*dt works well here. .2 also works.
			end
		else
			PreIndex = Index;
		end;

	end

	if j==iter
		disp('MBO did not converge'); 
	end
	[~,~,groups]=unique(Index);

	groups = greedy_cleanup(groups,W_sub,k,twom,gamma);

	Q_div= sum(sum(u.*(W_sub*u-(gamma/twom)*k*(k'*u))))/twom;
	Q_nodiv = (sum(sum(W_sub))-gamma*sum(k)^2/twom)/twom;
	dQ=Q_div-Q_nodiv;

end
