% This is a wrapper function for my_eigs.cpp. It "just works" in matlab.
% thread is needed to ensure that multiple instances can be run at once without trying to write the same file.
function [V,D] = my_eigs(W,gamma,k_ori,twom_ori,eigCount,thread,eig_type,opts)% no opts currently supported

	if nargin == 5
		eig_type = 'M';
		thread = 1;
	elseif nargin == 6
		eig_type = 'M';
	end

	if strcmp(eig_type,'CV')
		eig_type = 'M';
	elseif strcmp(eig_type,'MBO')
		eig_type = 'L';
	end

	N = size(W,2);

	if N > 500 && strcmp(eig_type,'M') %|| (strcmp(eig_type,'L') && N > 25) %this number was empirically determined for optimal speed. We exclude L for testing purposes. The 25 prevents pathologies, since Chris's solver was not designed to find all the eigenvalues of the matrix.

		% RC procedure

		write_W(W,gamma,k_ori,twom_ori,eigCount,thread,eig_type);
		olddir = cd('cpp_code');
		system(sprintf('./my_eigs.exe %d',thread));
		[V,D] = read_eigs(N,eigCount,thread);

		delete(sprintf('V_%d.txt',thread),...
			sprintf('D_%d.txt',thread),...
			sprintf('W_%d.dat',thread),...
			sprintf('params_%d.dat',thread),...
			sprintf('k_ori_%d.dat',thread));
		cd(olddir);

	elseif strcmp(eig_type,'L') && eigCount < N
		L = spdiags(sum(W,2),0,N,N) - W;
		[V,D] = eigs(L,eigCount,'sa');
	else

		% matlab's builtin

		k = sum(W,2);
		L = diag(k) - W;
		M = full(L+2*gamma/twom_ori*k_ori*k_ori');
		M = max(M,M');
		[V,D] = eig(M);
		V=V(:,1:eigCount);
		D=diag(D);
		D=D(1:eigCount);
	end
end
