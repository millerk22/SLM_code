% Based off of code writen by Hui Hu for the modularity MBO paper. Modified (a lot) by Zach Boyd.
% 
% W:  weight matrix
% gamma: resolution parameter for modularity
% Neig: number of eigenvectors/values of sparse Laplacian
% n: max number of clusters at each sub-partition step. (One can choose n=10)
% Output group: group assignments

%
function group = mainModMBOrecursionF(W,gamma,Neig,n_given, type,thread,verbose,randomInit,opts)

	if nargin < 6
		thread = 1;
	end
	if nargin < 7
		verbose = false;
	end

	tol=0.00000001; % tolerance for recursive sub-partition.
	run=1; % how many repeat runs for each recursion step.
	N=size(W,1);
	k=sum(W)';
	twom = full(sum(k));
	tosplit=true;
	group=ones(N,1);
	
	while any(tosplit)
		groupindex=unique(group);
		for ig=groupindex(tosplit(groupindex))'
			maxgroup=max(group);
			ind = find(group==ig);
			N_sub=length(ind);
			disp(N_sub);
			k_ori=k(ind);
			W_sub=W(ind,ind);

			%check for connectedness
			Index = find_components(W_sub);
			if nnz(Index) < N_sub && strcmp(type,'CV') % multiple connected components
					ic_sub = Index + 1;

					nhat = max(ic_sub); % check if modularity improves
					u = sparse(1:N_sub,ic_sub,ones(N_sub,1),N_sub,nhat);
					Q_div= sum(sum(u.*(W_sub*u-(gamma/twom)*k_ori*(k_ori'*u))))/twom;
					Q_nodiv = (sum(sum(W_sub))-gamma*sum(k_ori)^2/twom)/twom;
					dQ=Q_div-Q_nodiv;
			else %single connected component

				n=min([n_given round(sqrt(N_sub)) 11]); % n needs to be reasonable
				if(N_sub>1)
					n=max(n,2);
				end

				neig=min(2*n,N_sub);
				fast_cheat=false;
				if ~fast_cheat || N_sub<N
					[V,D] = my_eigs(W_sub,gamma,k_ori,twom,neig,thread,type);
					if N_sub==N
						save V_temp V D;
					end
				else
					disp('using fastcheat')
					load V_temp V D;
				end

				if randomInit
					[~,Index_1] = max(rand(N_sub,n),[],2);
					uo = sparse(1:N_sub,Index_1,1,N_sub,n);
					uo =full(uo);
				else
					disp('initializing with kmeans');
					uo=g2u(my_kmeans(V,n));
				end

				if strcmp(type, 'MBO')
					[ic_sub,dQ] = Modularity_MBO_sub(W_sub,uo,gamma,V,D,k_ori,twom);
				elseif strcmp(type,'CV')
					ssl=opts.ssl;
					if ssl.use_me==true
						ssl.use_me=false; disp('warning: ssl not currently supported in recursive mode');
					end
					if ssl.use_me
						disp('saving temp'); save temp ssl ind N_sub n; disp('save temp')
						ssl_int.bool = ssl.bool(ind);
						ssl_int.ind = find(ssl_int.bool);
						ssl_int.num = length(ssl_int.ind);
						ssl_int.g = ssl.g(ind);
						[~,~,temp] = unique(ssl_int.g(ssl_int.g>0));
						ssl_int.g(ssl_int.g>0) = temp;
						% A modification would need to be made here: if the subgraph contains more classes than we are 
						% partitioning into, then ssl_in.u will be bigger than we want. A simple solution is to select n
						% classes and use only the supervision from them, ignoring the rest. Ran out of time to implement
						% this.
						ssl_int.u = sparse(ssl_int.ind,ssl_int.g(ssl_int.g>0),1,N_sub,n); % Needs to be sparse for big nhat
						ssl_int.use_me=length(ssl_int.ind~=0);
						ssl_int.lam=ssl.lam;
					else
						ssl_int.use_me=false;
					end
					opts_int.ssl=ssl_int;
					opts_int.pseudospectral=opts.pseudospectral;
					[ic_sub,dQ] = Balanced_TV(W_sub,uo,gamma,V,D,k_ori,twom,opts_int);
				else
					warning('Unknown type in mainModMBOrecursion.m')
				end
			end

			if dQ >tol      
				for j=2:max(ic_sub)
					tempInd=(ic_sub==j);              
					group(ind(tempInd))= maxgroup+j-1;              
					tosplit(maxgroup+j-1)=(nnz(tempInd)>1);
				end
				tosplit(ig)=(nnz(ic_sub==1)>1);
			else
				tosplit(ig)=false;
			end
		end
	end

end
