function g = klnb(g,A,k,twom,gamma)

  if (size(g,1)>size(g,2)),
    g=g';
  end
  [~,~,g]=unique(g);
  u=g2u(g);

  N=length(g);
  nhat=max(g);
  if nhat==1
    return;
  end

  diagB = -gamma*k.^2/twom*ones(1,nhat);

  Nbrs=cell(N,1);
  for i=1:N
    Nbrs{i} = find(A(:,i));
  end

  pass_improvement=1;
  tol=10^(-6);
  while pass_improvement>tol %Outer loop until improvement not found inside
    2*pass_improvement/twom

    net_improvement=0;
    pass_improvement=0;

    indx=1:N;

    g_temp=g;
    u_temp = g2u(g_temp,nhat);

    Cut = A*u_temp;
    vol = k'*u_temp;

    no_improvement=0; % suggested in https://arxiv.org/pdf/0812.4073.pdf
    max_no_improvement=10*log2(nnz(A));

    while ~isempty(indx) %Move each node precisely once

      I = diagB... 
	+ Cut... 
	- sum(Cut.*u_temp,2)*ones(1,nhat)... 
	- gamma/twom*k*vol... 
	+ gamma/twom*k.*vol(g_temp)'*ones(1,nhat); 

      II=I(indx,:);%Restrict to available nodes
      II(u_temp(indx,:)==1)=-Inf;%Ignore current placement
      maxvalue=max(max(II)); %Identify best available switch

      [i,a]=find(II==maxvalue);
      ii=ceil(length(i)*rand(1)); %break ties randomly
      i=indx(i(ii)); a=a(ii);
      a_old = g_temp(i);

      indx(indx==i)=[]; %remove i from indx list of nodes to be moved

      g_temp(i)=a; %move node i to new group
      u_temp(i,a_old) = 0;
      u_temp(i,a)=1;

      Cut(Nbrs{i},a_old) = Cut(Nbrs{i},a_old) - A(Nbrs{i},i);
      Cut(Nbrs{i},a)     = Cut(Nbrs{i},a)     + A(Nbrs{i},i);

      vol(a_old) = vol(a_old) - k(i);
      vol(a)     = vol(a)     + k(i);

      net_improvement=net_improvement + maxvalue;

      if net_improvement > pass_improvement %Keep track of configuration with highest Q.
	g=g_temp;
	pass_improvement=net_improvement;
      else 
	no_improvement = no_improvement + 1;
	if no_improvement > max_no_improvement
	  indx=[]; % give up
	end
      end
    end
  end
end
