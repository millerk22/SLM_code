function write_W(W,gamma,k_ori,twom_ori, eigCount,thread,eig_type);

	W = W - diag(diag(W));

	[i,j,val] = find(W);
	W_dump = [i,j,val];

	% write everything
	fid = fopen(sprintf('cpp_code/W_%d.dat',thread),'w');

	if fid == -1 % this usually happens when we are in the wrong directory
		olddir = cd('..');
		fid = fopen(sprintf('cpp_code/W_%d.dat',thread),'w');
	end
	if fid ==-1 % Not sure what would cause this.
		cd(olddir); % Better cd back...
	end

	fprintf( fid,'%d %d %e\n', transpose(W_dump) );
	fclose(fid);

	% include important parameters
	fid = fopen(sprintf('cpp_code/params_%d.dat',thread),'w');
	N = size(W,2);
	twom_ori = full(twom_ori);
	n_nnz = nnz(W);
	fprintf( fid,'%d %d %e %e %d %s', n_nnz, N, gamma, twom_ori, eigCount, eig_type);
	fclose(fid);

	% write k_ori
	fid = fopen(sprintf('cpp_code/k_ori_%d.dat',thread),'w');
	k_ori = full(k_ori);
	fprintf( fid,'%d\n', k_ori');
	fclose(fid);
end
