function g = split(A,gamma,n_given,g,opts,to_split)
	global verbose;
	tol=0.00000001; % tolerance for recursive sub-partition.
	k=sum(A)';
	twom = full(sum(k));
	m = twom/2;
	u = g2u(g);
	Cut = u'*A*u;
	vol = k'*u;

	Q_old = Q_AIC(g,A,gamma,opts);

	to_split_ = (hist(g,max(g)) > 1);
	to_split = and(to_split,to_split_);

	disconnected = false; % a flag to see if disconnected communities are detected

	to_cont = true;
	while  to_cont == true;

		[~,ig] = max(vol.*to_split); % Try splitting the biggest ones first
		%ig = find(to_split,1); % Just pick the first one
		maxg=max(g);

		[ind,N_sub,A_sub] = get_subgraph(ig,g,A);
		if max(size(A_sub))==0; save temp; error('empty subgraph!'); end;
		if verbose; fprintf('split: trying to split a group of size: %d\n', N_sub); end

		% Propose a split
		ic_sub = get_components(A_sub);
		if max(ic_sub) > 1
			disconnected=true;
			if verbose; fprintf('found %d connected components\n',length(unique(ic_sub))); end;
		else

			n_split = get_n_split(n_given,N_sub);
			g0 	= get_g0(n_split,N_sub);

			ic_sub  = get_ic_sub(g0,A_sub,gamma,k(ind),m,opts.solver,opts); % Call the surface tension solver
			if verbose; fprintf('propose to split this group into %d communities\n',length(unique(ic_sub))); end;
		end

		g_temp = ic2g(g,ic_sub,ind);
		u_temp = g2u(g_temp);
		Cut_temp = get_Cut(g_temp,g,Cut,A);
		vol_temp = k'*u_temp;

		to_split_temp = get_to_split(to_split,ic_sub,ig,maxg);

		opts.to_split = to_split_temp;

		[g_temp,to_split_temp,Q_new,Cut_temp,vol_temp] = rc_merge(A,gamma,g_temp,Cut_temp,opts,k);

		%[g_temp,Cut_temp,vol_temp] = greedy_cleanup(g_temp,A,k,twom,gamma);
		Q_new=Qcut(g_temp,A,gamma,Cut_temp,vol_temp);

		if verbose; disp(sum(g2u(g_temp))); end;

		dQ = Q_new - Q_old;

		% Implement the split if it is good. Otherwise move on. Always split if the communities are not connected.
		if dQ < -tol || disconnected==true

			g = g_temp;
			to_split = to_split_temp;
			to_split(sum(g2u(g))<2) = 0; % Don't try to split singletons. (Allows later functions to make assumptions.)

			Q_old = Q_new;
			Cut = Cut_temp;
			vol = vol_temp;
			if verbose; fprintf('split: splitting proposal accepted\n'); end

		else
			if verbose; fprintf('split: splitting proposal rejected\n'); end
			to_split(ig)=false;
		end

		disconnected = false; % reset the flag
		to_cont = any(to_split); % see if more splitting is needed

		u = g2u(g);
		sz = sum(u(:,to_split));
		if verbose; disp(sz); end;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Auxilliary functions %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g_temp = ic2g(g,ic_sub,ind)
	g_temp = g;
	for j=2:max(ic_sub)
		g_temp(ind(ic_sub==j))= max(g)+j-1;              
	end
end

function [ind,N_sub,A_sub] = get_subgraph(ig,g,A);
	ind = find(g==ig);
	N_sub=length(ind);
	A_sub=A(ind,ind);
end

function n_split = get_n_split(n_given,N_sub);
	n_split=min(n_given,round(sqrt(N_sub))+1);
	if(N_sub>1)
		n_split=max(n_split,2);
	end
end

function g0 = get_g0(n_split,N_sub)
	g0 = randi(n_split,N_sub,1);
	g0(1:n_split) = 1:n_split;
end

function to_split = get_to_split(to_split,ic_sub,ig,maxg)
	for j=2:max(ic_sub)
		to_split(maxg+j-1)=(nnz(ic_sub==j)>1);
	end
	to_split(ig)=(nnz(ic_sub==1)>1);
end

function ic_sub = get_ic_sub(g0,A_sub,gamma,k,m,solver,opts)
	disp(size(A_sub,2));
	thread = opts.thread;
	type = 'M';
	eigCount = min(2*nnz(unique(g0)),length(g0));
	%disp('in debug mode');
	%disp(eigCount);
	%global zzz; % 94 may be optimal
	%disp(zzz)
	%eigCount=min(zzz,length(g0));
	if opts.pseudospectral
		fast_cheat=false;
		if ~fast_cheat || sum(sum(A_sub))<2*m || exist('V_temp.mat') ~= 2 %the second condition detects whether N_sub<N
			[V,D] = my_eigs(A_sub,gamma,k,2*m,eigCount,thread,type);
			if sum(sum(A_sub))==2*m
				save V_temp V D;
			end
		else
			disp('using fastcheat')
			load V_temp V D;
		end
	else
		V=0;
		D=[0 0];
	end
	[ic_sub,~] = Balanced_TV(A_sub,g2u(g0),gamma,V,D,k,2*m,opts);
	[~,~,ic_sub] = unique(ic_sub);
end

function Cut_temp = get_Cut(g_temp,g,Cut,A)

	cc = [unique(g(g~=g_temp)) (max(g)+1):max(g_temp)];
	oc = 1:max(g);
	oc(cc) = 0;
	oc = oc(oc ~= 0);

	cn = find(ismember(g,cc));
	on = find(not(ismember(g,cc)));

	u_temp = g2u(g_temp);

	Cut_temp(oc,oc) = Cut(oc,oc);
	Cut_temp(cc,cc) = u_temp(cn,cc)'*A(cn,cn)*u_temp(cn,cc);
	Cut_temp(cc,oc) = u_temp(cn,cc)'*A(cn,on)*u_temp(on,oc);
	Cut_temp(oc,cc) = Cut_temp(cc,oc)';
end
