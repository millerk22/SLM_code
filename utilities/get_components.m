% Takes a graph adjacency matrix and returns a vector indicating the connected component to which each node belongs. This is used to let you segment each component separately using the submatrix A( (tag == i), (tag==i) )
function tag = get_components(A)

	N = size(A,2);
	tag = zeros(N,1);
	count = 1;
	while nnz(tag==0)>0
		seed = find(tag==0,1);
		u = get_one_component(A,seed);
		tag(u==1) = count;
		count = count + 1;
	end
		
