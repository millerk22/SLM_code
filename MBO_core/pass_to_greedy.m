function g = pass_to_greedy(g,A,k,twom,gamma,thread)

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

	write_vec(g,'g',thread);
	write_spr(A,'A',thread);
	write_vec(k,'k',thread);
	write_dbl(twom,'twom',thread);
	write_dbl(gamma,'gamma',thread);
	write_vec(diagB,'diagB',thread);
	mat(Cut,'Cut',thread);
	write_vec(vol,'Cut',thread);

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

	g = read_vec('g',thread);

	disp('end cleanup');
end

function V = read_vec(name,thread)

	fileID =fopen(sprintf('%s_%d.txt',name,thread));

	if fileID == -1
		error('file.exe failed to open');
	end

	formatspec = '%f';
	sizeA = [1 Inf];
	V = fscanf(fileID,formatspec,sizeA);
	fclose(fileID);

end

function write_vec(V,name,thread)

	fid = fopen(sprintf('%s_%d.dat',name,thread),'w');
	V=full(V);
	fprintf( fid,'%d\n', V');
	fclose(fid);

end

function write_dbl(V,name,thread)
	fid = fopen(sprintf('%s_%d.dat',name,thread),'w');
	fprintf( fid,'%d', V);
	fclose(fid);
end

function write_spr(A,name,thread)
	fid = fopen(sprintf('%s_%d.dat',name,thread),'w');

	fprintf(fid,'%d %d %f\n', transpose(A) );
	fclose(fid);
end

function write_mat(A,name,thread)
	fname = fopen(sprintf('%s_%d.dat',name,thread),'w');
	dlmwrite(fname,A);
end
