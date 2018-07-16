function params = NEW_get_data(str,opts);

  % defaults
  nhat = 2;
  eigCount = 2*nhat;
  gamma = 1;
 

  switch str
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'pnnl_test'
      load W_test.mat
      tag=W(:,1);
      W=W(:,2:size(W,2));
      W = max(W,W');
      W = W - diag(diag(W));
      N=size(W,2);%Number of Nodes in the graph
      gamma=1;
      nhat=2;% number of communties
      eigCount=nhat*2;
      have_nhat = true;
      have_ssl=true;
      
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
      error('Incorrect problem name.');

  end
% 
%   Standard postprocessing
 
  fulltag=tag;
  params.fulltag=fulltag;
  [~,~,tag] = unique(tag);
%   [W,tag,ind] = get_giant_component(W,tag);
  ind= 1: size(W,1);
  params.A = W;
  params.eigCount = eigCount;
  params.nhat = nhat;
  params.gamma = gamma;
  params.tag = tag;
  params.have_nhat = have_nhat;
  params.have_ssl=have_ssl;
  params.sstag = tag;
  params.ind=ind;
  params.fulltag=fulltag;
