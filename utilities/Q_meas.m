%Calculates modularity of a network partition
%Q is the modularity
%group is a vector containing th group assignment of each node
%W is the graph adjacency matrix
function Q = Q_meas(group, W, gamma)

	[~,~,group] = unique(group);%makes sure the numbers are as expected.

	k = sum(W)';
	twom = sum(k);

	Q = 0;
	for i=1:max(group)
		iind=find(group==i);
		Q = Q + sum(sum(W(iind,iind)))-gamma*sum(k(iind))^2/twom;
	end
	Q = Q/twom;

	Q = full(Q);

end
