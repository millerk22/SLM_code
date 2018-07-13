%Reads in the spectral data created through c++ routines.
function [V,D] = read_eigs(N,Neigs,thread)

	V=zeros(N,Neigs);
	D=zeros(Neigs,1);

	fileID =fopen(sprintf('V_%d.txt',thread));

	if fileID == -1
		cd('..');
		error('my_eigs.exe failed to produce eigenvectors');
	end

	formatspec = '%e';
	sizeA = [1 Inf];
	V = fscanf(fileID,formatspec,sizeA);
	fclose(fileID);

	V = reshape(V,[N,Neigs]);

	fileID =fopen(sprintf('D_%d.txt',thread));
	formatspec = '%e';
	sizeA = [1 Inf];
	D = fscanf(fileID,formatspec,sizeA);
	D=D';
	fclose(fileID);
end
