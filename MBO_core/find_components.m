% This function takes a graph adjacency matrix in which it is known that there are at least two connected components and returns on connected component. The approach is to repeated walks until a steady state is reached.
% We return a vector of 0s and 1s indicating the separation.

function u = find_components(W)
	N = size(W,1);
	u = zeros(N,1);
	u_old=u;
	u(1) = 1;
	while ~isequal(u,u_old)
		u_old = u;
		u = u + W*u; % do a random walk
		u = (u>0); % enforce binary
	end
end
