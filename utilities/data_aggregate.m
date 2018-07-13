imax=20;
addpath(genpath(pwd));

types={'louvain','CVR'};
tmax=length(types);

for d=1:2
  dt=50*d;
  dt
  params=get_data(sprintf('pnnl_%d',dt));
  N=length(params.fulltag);
  QQ=zeros(tmax,imax);
  gg=zeros(tmax,imax,N);
  for t=1:tmax
    type=types{t};
    type
    for i=1:imax
      i
      fname=sprintf('pnnl_%d_%s_%d.mat',dt,type,i);
      if exist(fname,'file')==2
	load(fname);
	convert_to_old_indices;

	QQ(t,i)=Qg;
	gg(t,i,:)=g;
    end
  end
  fname=sprintf('pnnl_%d_%s_group.csv',dt,type);
  csvwrite(fname,squeeze(gg(t,:,:)));
  fname=sprintf('pnnl_%d_%s_modularity.csv',dt,type);
  csvwrite(fname,QQ(t,:));
end
end
