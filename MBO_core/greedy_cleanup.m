% CCut is the between-class cuts
% Cut  is the between-node  cuts
function [g,CCut,vol] = greedy_cleanup(g,A,k,twom,gamma,CCut)

	disp('cleanup');

	if (size(g,1)>size(g,2)),
		g=g';
	end
	[~,~,g]=unique(g);
	u=g2u(g);

	N=length(g);
	nhat=max(g);

	diagB = -gamma*k.^2/twom*ones(1,nhat);

	Nbrs=cell(N,1);
	for i=1:N
		Nbrs{i} = find(A(:,i));
	end

	Cut = A*u;
	vol = k'*u;

	tol=10^(-4);
	pass_improvement=Inf;
	while 2*pass_improvement/twom>tol %Outer loop until improvement not found inside

		pass_improvement=0;

		for i=1:N

			a_old=g(i);

			I = Cut(i,:)-gamma/twom*k(i)*vol;

			[maxvalue,a]=max(I);

			if a ~= a_old

				g(i)=a;

				pass_improvement=pass_improvement + maxvalue + diagB(i) - Cut(i,a_old) + gamma/twom*k(i)*vol(a_old);

				Cut(Nbrs{i},a_old) = Cut(Nbrs{i},a_old) - A(Nbrs{i},i);
				Cut(Nbrs{i},a)     = Cut(Nbrs{i},a)     + A(Nbrs{i},i);

				vol(a_old) = vol(a_old) - k(i);
				vol(a)     = vol(a)     + k(i);
			end
		end

		2*pass_improvement/twom
	end

	CCut = g2u(g,nhat)'*Cut;

	disp('end cleanup');
end
