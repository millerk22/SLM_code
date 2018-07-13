function group = my_kmeans(V,nhat) %Standard kmeans

	group=-1;
	count = 1;
	while group(1)==-1 && count < 4
		group=my_kmeans_int(V,nhat);
		count=count+1;
	end

	if count==4;
		warning('kmeans failed three times. continuing with random guess.')

		N=size(V,1);
		group = randi(nhat,1,N);
		x=min(nhat,N);
		group(1:x)=1:x;
	end

end

function group = my_kmeans_int(V,nhat) %Standard kmeans
	V=V(:,1:nhat); % in case more are provided

	dim = size(V,2);
	mu = zeros(nhat,dim);

	N = size(V,1);

	%Random initial assignments
	group = randi(nhat,N,1);
	group(1:nhat) = 1:nhat; % Each group must be non-empty

	tol = max(0.05,4/N); % smaller graphs need a higher tolerance to avoid cycling
	sim = 0;

	while sim < 1-tol

		sim
		old_group = group;

		% Reassign the centers
		for i=1:nhat
			if nnz(group==i)>0
				mu(i,:) = mean(V(group==i,:));
			else
				disp('kmeans failed once. retrying.')
				g=-1;
				return;
			end
		end

		for j=1:nhat
			temp(:,j) = sum((V - ones(N,1)*mu(j,:)).^2,2);
		end
		[~,group] = min(temp,[],2);

		sim = purityMeas(group,old_group);
		if exist('old_old_group')
			sim = max(sim,purityMeas(group,old_old_group));
		end
		old_old_group=old_group;

	end
end
