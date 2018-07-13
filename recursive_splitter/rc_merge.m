% Q_opts can be empty if needed.
function [g,to_split,Q,Cut,vol] = rc_merge(A,gamma,g,Cut,Q_opts,k)

	global verbose;
	if verbose; disp('entering rc_merge'); end;
	if verbose; fprintf('%d communities found\n',max(g)); end;
	if verbose; disp(sum(g2u(g))); end;

	to_split = [];
	if isfield(Q_opts,'to_split') 
		to_split = Q_opts.to_split; 
		to_split = find(to_split); % The split and merge codes use a different format for this part
	end;

	twom = full(sum(k));

	merge_found = true;
	u = g2u(g);
	vol = k'*u;
	Q_opt = Q_AIC(g,A,gamma,Q_opts,Cut,vol);
	merge_list = Cut - diag(diag(Cut));
	while merge_found

		Q_old = Q_opt; % Value found in previous iteration

		opt_improvement = 0;
		merge_found = false;

		while nnz(merge_list) > 0
			[i,~] = get_pair(merge_list);
			while nnz(merge_list(i,:)) > 0

				j = find(merge_list(i,:),1);
				i_old = i;
				if j<i %While strictly unnecessary, this gives easier comparison with output of Qmerge, since it does this internally 
					t = i;i = j;j = t;
				end
				%[i,j] = get_pair(merge_list);

				%%%%%%%%%%% Get improvement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				g_temp = g;
				g_temp(g_temp == j) = i;
				g_temp(g_temp>j) = g_temp(g_temp>j) - 1; % Could use unique here, but I think it is slower

				[Cut_temp,vol_temp,Q_new] = Qmerge(Cut,vol,i,j,twom,gamma,Q_opts);
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

				temp_improvement = Q_old - Q_new;

				tol = 0.000001;
				if temp_improvement > opt_improvement + tol;

					i_opt = i;
					j_opt = j;
					g_opt = g_temp;
					Cut_opt = Cut_temp;
					vol_opt = vol_temp;
					Q_opt = Q_new;

					opt_improvement = temp_improvement;

					merge_found = true;
				end

				merge_list(i,j) = 0;
				merge_list(j,i) = 0;

				i = i_old; % Potentially undo the swap earlier
			end
		end

		if merge_found
			ratio = vol(i)/vol(j);
			[g,to_split,merge_list,Cut,vol] = do_merge(i_opt,j_opt,to_split,g_opt,Cut_opt,vol_opt,Q_opt,ratio);
		end

		Q = Q_old;

	end

	% Change format to be compatible with split
	temp = zeros(max(g),1);
	temp(to_split) = 1;
	to_split = (temp==1)';
end

function [i,j] = get_pair(merge_list)
	global verbose;
	[i,j] = find(merge_list == max(max(merge_list)),1);
	if j<i %strictly unnecessary, makes for easier comparison with output of Qmerge, since it does this internally
		t = i;i = j;j = t;
	end
end

function [g,to_split,merge_list,Cut,vol] = do_merge(i,j,to_split,g_temp,Cut_temp,vol_temp,Q_new,ratio)
	global verbose;

	if verbose; fprintf('\nrc_merge: merging %d and %d\n',i,j); end

	g = g_temp;

	Cut = Cut_temp;
	vol = vol_temp;

	merge_list = Cut - diag(diag(Cut));

	count =1;
	Q_old = Q_new;

	to_split = to_split(to_split ~= j);
	if ratio < 0.9 && ratio > 0.1 % only resplit i if the volume of i and j are comparable
		to_split = [to_split i];
	end
	to_split(to_split > j) = to_split(to_split>j) - 1;
	%to_split = [to_split find(Cut(i - (i>j),:))];
	to_split = unique(to_split);
end
