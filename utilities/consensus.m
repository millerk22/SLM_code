for iter = 2:nhat
	refgroup = groups(:,1);
	oldgroup = groups(:,iter);
	newgroup = groups(:,iter);

	for al = 1:nhat
		
		inds = find(oldgroup == al);
		newind = mode(refgroup(inds));
		%newind = newind(1);
		newgroup(inds) = newind;

	end

	groups(:,iter) = newgroup;
end

consensusgroup = mode(groups,2);
