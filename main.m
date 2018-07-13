if ~exist('benchmark_flag') || ~benchmark_flag
	problem = 'karate'
	type = 'CV' % CVR, CV, louvain, newman, clauset
	already_have_eigs = false % in case we computed them on a previous run
	deterministic = false
	thread = 22430544;
	semisupervised=false
	pseudospectral=true
	recurse=false
	randomInit=true
end

addpath(genpath(pwd));

if deterministic == true
	seed=81;
	my_rng(seed);
	opts.seed=seed;
else
	seed=randi(100000,1);
	opts.dummy = true;
end
fid = fopen(sprintf('cpp_code/seed_%d.dat',thread),'w');
fprintf( fid,'%d',seed);
fclose(fid);
opts.thread = thread;

params = get_data(problem,opts);

A = params.A;
eigCount=params.eigCount; %Use 20, but fewer also works.
gamma=params.gamma;
tag = params.tag;
have_nhat=params.have_nhat;
if have_nhat && recurse
	disp('using recursion even though this dataset has a known number of communities');
end
recurse = (~have_nhat || recurse || strcmp(type,'CVR'));
nhat=params.nhat;
N = length(tag);
if ~params.have_ssl && semisupervised
	disp('semisupervision requested but is unavailable for this dataset');
end
semisupervised = params.have_ssl && semisupervised; % both the data and the caller must approve;

ssl.use_me=semisupervised;
if semisupervised
	disp('semisupervised');
	frac = 0.1;
	ssl.bool = rand(1,N) > 1-frac;
	ssl.bool = ssl.bool .* (params.sstag' > 0);
	ssl.ind = find(ssl.bool);
	ssl.num = length(ssl.ind);
	ssl.g = params.sstag.*ssl.bool';
	ssl.u = sparse(ssl.ind,params.sstag(ssl.ind),1,N,nhat); % Needs to be sparse in case nhat is big
	ssl.lam = 10;
end

N = size(A,2);
L=spdiags(sum(A)',0,N,N)-A;%unnormalized laplacian
k=full(sum(A,2));
twom=sum(k);

if (~recurse || ~randomInit) && ~already_have_eigs

	if ~pseudospectral 
		eigCount = nhat;
	end
	disp('getting eigs')
	tic
	if strcmp(type,'CV') || strcmp(type,'CVR')
		[V,D] = my_eigs(A,gamma,k,twom,eigCount,thread,type);
		D=diag(D);
	elseif strcmp(type,'MBO')
		[V,D] = my_eigs(A,gamma,k,twom,eigCount,thread,type);
		D=diag(D);
	end
	eig_time = toc;
end

if randomInit && ~recurse
	[~,Index] = max(rand(N,nhat),[],2);
	uo = full(sparse(1:N,Index,ones(N,1),N,nhat));
elseif ~recurse
	disp('initializing with kmeans')
	uo=g2u(my_kmeans(V,nhat));
end

%% Modularity_MBO minimization
tic
disp('calling main algorithm');
opts.ssl=ssl;
opts.pseudospectral=pseudospectral;
if strcmp(type,'MBO') && ~recurse
	g = Modularity_MBO(A,uo,gamma,V,D);
elseif strcmp(type,'CV') && ~recurse
	disp('no recursion')
	g = Balanced_TV(A,uo,gamma,V,D,k,twom,opts);
elseif strcmp(type,'CVR')
	opts.nhat=nhat;
	opts.lam=0.0; disp('setting lam=0')
	opts.solver = type;
	[g,~] = recursive_splitter(A,gamma,min(nhat,4),opts);
elseif  (strcmp(type,'MBO') || strcmp(type,'CV'))
	disp('using recursion')
	g = mainModMBOrecursionF(A,gamma,eigCount,nhat,type,thread,true,randomInit,opts);
elseif strcmp(type,'louvain') || strcmp(type,'newman') || strcmp(type,'clauset')
	if ~exist('thread')
		thread = 1;
	end
	g = igraph_call(type,A,gamma,thread);
else
	warning('Invalid type in main');
end
times = toc;

%Check performance
Qtag = Q_meas(tag,A,gamma)
Qg = Q_meas(g,A,gamma)
pur= purityMeas(g,tag)
nmi = nmiMeas(g,tag)
