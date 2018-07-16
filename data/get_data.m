function params = get_data(str,opts);

  % defaults
  nhat = 4;
  eigCount = 2*nhat;
  gamma = 1;
  have_nhat = false;
  have_ssl=false; % can "tag" be used as ground truth for semi-supervised learning?

  switch str
    case 'karate'
      [W,tag] = get_karate();
      N=34;
      nhat=2;
      have_nhat = true;
      have_ssl=true;
    case 'twoMoons'
      load modsamples2000
      eigCount=20; %Use 20, but fewer also works.
      N=2000;
      maxIter=35; %Use up to 35
      gamma=.2;
      nhat=2;
      %iter = randi(35,1);
      iter = 5;
      W=Wm((iter-1)*N+1:iter*N,:);
      tag = squeeze(tagm(iter,:))';
      have_nhat = true;
      have_ssl=true;
    case 'MNIST'

      gamma=0.5;
      nhat=11;
      eigCount = 100; % number of eigenvectors/values of sparse Laplacian. 100 works well.
      load('weight_mnist70k')
      tag = Label;
      W=w;
      have_nhat = true;
      have_ssl=true;

      %		case 'LFR'
      %
      %			eigCount=40;%Number of eigenvectors
      %			N = 50000;
      %			gamma = 15; %does give about the right number of communities
      %			nhat = 40;
      %			mix = 1;
      %			%iter = 1;
      %			iter = randi(4,1);
      %			[W,tag]=getdata50(mix,iter);

    case 'LFR'

      opts.N = 50000;opts.k_av = 20;opts.k_max = 50; opts.lam = 0.001;
      N = opts.N;
      %make_data(opts);
      %[W,tag]=lfr_getdata(N);
      load W_lfr_50000;
      [W,tag] = get_giant_component(W,tag);
      N = size(W,2);
      nhat = max(tag);
      gamma = 15;
      eigCount = 20;
      have_nhat = false;

      have_ssl=true;
    case 'urban'
      load data;
      tic;
      %W=kdtreeweight(H,10); % No NL means
      %W=kdtreeweight_3(H,10); %3x3 window
      load W_urban;
      W=(W+W')/2;
      toc

      N = size(W,1);

      %resolution parameter
      %gamma = 1; % Produces many classes. Generally appealing.
      gamma = .1; % Produces 5 classes. Generally lower quality.
      nhat=6; % Works well for recursion
      eigCount=2*nhat;%Number of eigenvectors

      W = W-diag(diag(W));
      tag = ones(N,1);
      have_nhat = true;
    case 'cuprite'
      load W_cuprite;
      % from http://www.escience.cn/people/feiyunZHU/Dataset_GT.html

      N = size(W,1);

      %resolution parameter
      %gamma = 1; % Produces many classes. Generally appealing.
      gamma = .1; % Produces 5 classes. Generally lower quality.
      nhat=12; 
      eigCount=2*nhat;%Number of eigenvectors

      W = W-diag(diag(W));
      tag = ones(N,1);
      have_nhat = true;

    case 'jasperRidge'
      load W_jasperRidge;
      % from http://www.escience.cn/people/feiyunZHU/Dataset_GT.html

      N = size(W,1);

      %resolution parameter
      %gamma = 1; % Produces many classes. Generally appealing.
      gamma = .1; % Produces 5 classes. Generally lower quality.
      nhat=4; % Works well for recursion
      eigCount=2*nhat;%Number of eigenvectors

      W = W-diag(diag(W));
      tag = ones(N,1);
      have_nhat = true;


    case 'samson_1'
      load W_samson_1;
      % from http://www.escience.cn/people/feiyunZHU/Dataset_GT.html

      N = size(W,1);

      %resolution parameter
      %gamma = 1; % Produces many classes. Generally appealing.
      gamma = .1; % Produces 5 classes. Generally lower quality.
      nhat=3; % Works well for recursion
      eigCount=2*nhat;%Number of eigenvectors

      W = W-diag(diag(W));
      tag = ones(N,1);
      have_nhat = true;
    case 'flc1'
      load W_flc1;
      %https://engineering.purdue.edu/~biehl/MultiSpec/hyperspectral.html

      N = size(W,1);

      %resolution parameter
      %gamma = 1; % Produces many classes. Generally appealing.
      gamma = .1; % Produces 5 classes. Generally lower quality.
      nhat=3; % Works well for recursion
      eigCount=2*nhat;%Number of eigenvectors

      W = W-diag(diag(W));
      tag = ones(N,1);
      have_nhat = true;
    case 'paviaU'
      load W_paviaU;

      N = size(W,1);

      %resolution parameter
      %gamma = 1; % Produces many classes. Generally appealing.
      gamma = .1; % Produces 5 classes. Generally lower quality.
      nhat=9; % Works well for recursion
      eigCount=2*nhat;%Number of eigenvectors

      W = W-diag(diag(W));
      tag = ones(N,1);
      have_nhat = true;
    case 'salinas'
      load W_salinas;

      N = size(W,1);

      %resolution parameter
      %gamma = 1; % Produces many classes. Generally appealing.
      gamma = .1; % Produces 5 classes. Generally lower quality.
      nhat=6; % Works well for recursion
      eigCount=2*nhat;%Number of eigenvectors

      W = W-diag(diag(W));
      tag = ones(N,1);
      have_nhat = true;
    case 'salinas_1'
      load W_salinas_1;

      N = size(W,1);

      %resolution parameter
      %gamma = 1; % Produces many classes. Generally appealing.
      gamma = .1; % Produces 5 classes. Generally lower quality.
      nhat=16; % Works well for recursion
      eigCount=2*nhat;%Number of eigenvectors

      W = W-diag(diag(W));
      tag = ones(N,1);
      have_nhat = true;
    case 'plumes7'

      load W_7_3;
      N = size(W,1);
      gamma = 1; % 0.5 works better for some recursive work. use 1 for compatibility with igraph
      nhat = 6;
      Neig = 2*nhat;
      tag = ones(N,1);
      have_nhat = true;

    case 'plumes40'

      load W_40_3 W; %loads W, the adjacency matrix. For memory reasons, this was easiest to compute on Joshua.
      N = size(W,1);
      gamma = 1; % 0.5 works better for some recursive work. use 1 for compatibility with igraph
      nhat = 6;
      Neig = 2*nhat;
      tag = ones(N,1);
      have_nhat = true;

    case 'Caltech'
      load Caltech36; 
      W=A;
      nhat=8; % W, expected to have the community structure of houses

      N = size(W,2);
      tag = local_info(:,5);
      [~,~,tag] = unique(tag);
      [W,tag,ind] = get_giant_component(W,tag);
      N = size(W,2);
      have_nhat = true;
      have_ssl=true;
    case 'Princeton'
      load Princeton12; 
      W=A;
      nhat=4; % W, expected to have the community structure of years

      N = size(W,2);
      tag = local_info(:,6);
      [~,~,tag] = unique(tag);
      [W,tag,ind] = get_giant_component(W,tag);
      N = size(W,2);
      have_nhat = true;
      have_ssl=true;
    case 'Penn'
      load Penn94; 
      W=A;
      nhat=8; % W

      N = size(W,2);
      tag = local_info(:,6);
      [~,~,tag] = unique(tag);
      [W,tag,ind] = get_giant_component(W,tag);
      N = size(W,2);
      have_nhat = true;
      have_ssl=true;
    case 'Harvard'
      load Harvard1
      fb_script;
    case 'Columbia'
      load Columbia2
      fb_script;
    case 'Stanford'
      load Stanford3
      fb_script;
    case 'Yale'
      load Yale4
      fb_script;
    case 'Cornell'
      load Cornell5
      fb_script;
    case 'Dartmouth'
      load Dartmouth6
      fb_script;
    case 'UPenn'
      load UPenn7
      fb_script;
    case 'MIT'
      load MIT8
      fb_script;
    case 'NYU'
      load NYU9
      fb_script;
    case 'BU'
      load BU10
      fb_script;
    case 'Brown'
      load Brown11
      fb_script;
    case 'Berkeley'
      load Berkeley13
      fb_script;
    case 'Duke'
      load Duke14
      fb_script;
    case 'Georgetown'
      load Georgetown15
      fb_script;
    case 'UVA'
      load UVA16
      fb_script;
    case 'BC'
      load BC17
      fb_script;
    case 'Tufts'
      load Tufts18
      fb_script;
    case 'Northeastern'
      load Northeastern19
      fb_script;
    case 'UIllinois'
      load UIllinois20
      fb_script;
    case 'UF'
      load UF21
      fb_script;
    case 'Wellesley'
      load Wellesley22
      fb_script;
    case 'Michigan'
      load Michigan23
      fb_script;
    case 'MSU'
      load MSU24
      fb_script;
    case 'Northwestern'
      load Northwestern25
      fb_script;
    case 'UCLA'
      load UCLA26
      fb_script;
    case 'Emory'
      load Emory27
      fb_script;
    case 'UNC'
      load UNC28
      fb_script;
    case 'Tulane'
      load Tulane29
      fb_script;
    case 'UChicago'
      load UChicago30
      fb_script;
    case 'Rice'
      load Rice31
      fb_script;
    case 'WashU'
      load WashU32
      fb_script;
    case 'UC33'
      load UC33
      fb_script;
    case 'UCSD'
      load UCSD34
      fb_script;
    case 'USC'
      load USC35
      fb_script;
    case 'UCSB'
      load UCSB37
      fb_script;
    case 'Rochester'
      load Rochester38
      fb_script;
    case 'Bucknell'
      load Bucknell39
      fb_script;
    case 'Williams'
      load Williams40
      fb_script;
    case 'Amherst'
      load Amherst41
      fb_script;
    case 'Swarthmore'
      load Swarthmore42
      fb_script;
    case 'Wesleyan'
      load Wesleyan43
      fb_script;
    case 'Oberlin'
      load Oberlin44
      fb_script;
    case 'Middlebury'
      load Middlebury45
      fb_script;
    case 'Hamilton'
      load Hamilton46
      fb_script;
    case 'Bowdoin'
      load Bowdoin47
      fb_script;
    case 'Vanderbilt'
      load Vanderbilt48
      fb_script;
    case 'Carnegie'
      load Carnegie49
      fb_script;
    case 'UGA'
      load UGA50
      fb_script;
    case 'USF'
      load USF51
      fb_script;
    case 'UCF'
      load UCF52
      fb_script;
    case 'FSU'
      load FSU53
      fb_script;
    case 'GWU'
      load GWU54
      fb_script;
    case 'Johns Hopkins'
      load Johns55
      fb_script;
    case 'Syracuse'
      load Syracuse56
      fb_script;
    case 'Notre Dame'
      load 'Notre Dame57.mat'
      fb_script;
    case 'Maryland'
      load Maryland58
      fb_script;
    case 'Maine'
      load Maine59
      fb_script;
    case 'Smith'
      load Smith60
      fb_script;
    case 'UC64'
      load UC64
      fb_script;
    case 'VIllanova'
      load VIllanova
      fb_script;
    case 'Virginia'
      load Virginia63
      fb_script;
    case 'UC61'
      load UC61
      fb_script;
    case 'Cal'
      load Cal65
      fb_script;
    case 'Mississippi'
      load Mississippi66
      fb_script;
    case 'Mich'
      load Mich67
      fb_script;
    case 'UCSC'
      load UCSC68
      fb_script;
    case 'Indiana'
      load Indiana69
      fb_script;
    case 'Vermont'
      load Vermont70
      fb_script;
    case 'Auburn'
      load Auburn71
      fb_script;
    case 'USFCA'
      load USFCA72
      fb_script;
    case 'Wake'
      load Wake73
      fb_script;
    case 'Santa'
      load Santa74
      fb_script;
    case 'American'
      load American75
      fb_script;
    case 'Haverford'
      load Haverford76
      fb_script;
    case 'William'
      load William77
      fb_script;
    case 'MU'
      load MU78
      fb_script;
    case 'JMU'
      load JMU79
      fb_script;
    case 'Texas'
      load Texas8080
      fb_script;
    case 'Simmons'
      load Simmons81
      fb_script;
    case 'Bingham'
      load Bingham82
      fb_script;
    case 'Temple'
      load Temple83
      fb_script;
    case 'Texas84'
      load Texas84
      fb_script;
    case 'Vassar'
      load Vassar85
      fb_script;
    case 'Pepperdine'
      load Pepperdine86
      fb_script;
    case 'Wisconsin'
      load Wisconsin87
      fb_script;
    case 'Colgate'
      load Colgate88
      fb_script;
    case 'Rutgers'
      load Rutgers89
      fb_script;
    case 'Howard'
      load Howard90
      fb_script;
    case 'UConn'
      load UConn91
      fb_script;
    case 'UMass'
      load UMass92
      fb_script;
    case 'Baylor'
      load Baylor93
      fb_script;
    case 'Tennessee'
      load Tennessee95
      fb_script;
    case 'Lehigh'
      load Lehigh96
      fb_script;
    case 'Oklahoma'
      load Oklahoma97
      fb_script;
    case 'Reed'
      load Reed98
      fb_script;
    case 'Brandeis'
      load Brandeis99
      fb_script;
    case 'Trinity'
      load Trinity100
      fb_script;
    case 'Stelzl'
      W = bio_read(str);
      tag = ones(size(W,1),1);
    case 'YeastL'
      W = bio_read(str);
      tag = ones(size(W,1),1);
    case 'figeys'
      W = bio_read(str);
      tag = ones(size(W,1),1);
    case 'moreno'
      W = bio_read(str);
      tag = ones(size(W,1),1);
    case 'pdz'
      W = bio_read(str);
      tag = ones(size(W,1),1);
    case 'reactome'
      W = bio_read(str);
      tag = ones(size(W,1),1);
    case 'vidal'
      W = bio_read(str);
      tag = ones(size(W,1),1);
    case 'pnnl_50'
      load W_50.mat
      W=A;
      tag=ones(size(W,1),1);
    case 'pnnl_si'
      load W_pnnlv4_sig.mat
      tag=ones(size(W,1),1); 
    case 'pnnl_v4'
      load W_v4.mat
      tag=ones(size(W,1),1);
    case 'pnnl_fvb4'
      load W_b44fv.mat
      tag=ones(size(W,1),1);
    case 'pnnl_v4_B0'
      load W_B0.mat
      tag=ones(size(W,1),1);
    case 'pnnl_100'
      load W_100.mat
      W=A;
      tag=ones(size(W,1),1);
    case 'sig1'
      load W_S1.mat 
      tag=ones(size(W,1),1);
    case 'pnnl_test'
      load W_test2.mat
      tag=ones(size(W,1),1);
    case 'sig2'
      load W_S2.mat
      tag=ones(size(W,1),1);	
    case 'opioid_eucli'
      load A_eudli
      W=A;
      tag=ones(size(W,1),1);
    case 'opioid_wm'
      load A_wmd
      W=A;
      tag=ones(size(W,1),1);
    case '50KComposite'
      [W,tag] = get_50Kdata_2();
      N = 109084;
      %nhat = 100;
      have_nhat = true;
      have_ssl = true;
    case '100KComposite'
      [W,tag] = get_100Kdata();
      N = 218565;
      have_nhat = true;
      have_ssl = true;
    case 'PNNL4-B0'
      [W,tag,N] = get_PNNL_B0();
      nhat = 2;        % we want to have just 2 communities, signal nodes and non-signal nodes
      have_nhat = true; 
      have_ssl = true; % hold on, we haven't implemented this with ssl yet I thought...?
    otherwise
      error('Incorrect problem name.');

  end

  % Standard postprocessing
  W = max(W,W');        % make symmetric matrix
  W = W - diag(diag(W));    % ensure that the diagonal is 0 (no self loops)
  fulltag=tag;          % the fulltag is the full correct labeling of node classes
  [~,~,tag] = unique(tag); 
  [W,tag,ind] = get_giant_component(W,tag); % get largest connected component 
                            % and identify the labelings (tag) and indices (ind) in the
                            % original graph.

  % define the parameters for our problems
  params.A = W;
  params.eigCount = eigCount;
  params.nhat = nhat;
  params.gamma = gamma;
  params.tag = tag;
  params.have_nhat = have_nhat;
  params.have_ssl=have_ssl;
  if exist('sstag')
    params.sstag = sstag;
  else
    params.sstag = tag;
  end
  params.ind=ind;
  params.fulltag=fulltag;
