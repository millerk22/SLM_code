% Takes an adjacency matrix and returns the connected component containing the node given by seed
% We return a vector of 0s and 1s indicating the separation.

function u = get_one_component(A,seed)

	if nargin == 1
		seed = 1;
	end

	N = size(A,1);
	u = zeros(N,1);
	u_old=u;
	u(seed) = 1;
	while ~isequal(u,u_old)
		u_old = u;
		u = u + A*u; % do a random walk
		u = (u>0); % enforce binary
	end
end
