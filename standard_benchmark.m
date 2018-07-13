% Runs each method 3 times on each problem

benchmark_flag = true;
deterministic=false; % this may not do what you think it does
sb_iter_max=40;
global verbose
verbose=false;
semisupervised=false; %default false
pseudospectral=true; %default false
recurse=false;%default false
randomInit=true;%default false

if ~exist('list_flag') || ~list_flag % Am I being called from another script
	problem_list = {'Caltech'};
	solver_list = {'CVR'};
	thread = 8035234;
end

nprobs = length(problem_list);
nsolve = length(solver_list);

result_names = {'Q','T','Te','NMI','PUR'}; % names used by this script
local_names = {'Qg','times','eig_time','nmi','pur'}; % names as used by main
ops = {@max,@median,@median,@max,@max};
nresults = length(result_names);
for i=1:nresults
	results{i} = zeros(nprobs,nsolve);
end

pcount = 0;
for problem_singleton = problem_list
	pcount = pcount + 1;

	scount = 0;
	for solver = solver_list
		scount = scount + 1;

		for i=1:nresults
			sb_results{i} = zeros(sb_iter_max,1);
		end

		already_have_eigs = false;
		for sb_iter = 1:sb_iter_max
			problem = problem_singleton{1};
			type = solver{1};
			fprintf('\nSolving %s with %s, iteration %d\n',problem,type,sb_iter)
			if strcmp(type,'CV') || strcmp(type,'MBO') || strcmp(type,'CVR') || (sb_iter==1) % Some methods are deterministic
				main % uses already_have_eigs
				%fname = sprintf('%s_%s_%d.mat',problem,type,sb_iter);
				%save(fname,'g','Qg');
			end

			for i=1:nresults
				if exist(local_names{i})
					sb_results{i}(sb_iter) = eval(local_names{i});
				end
			end
			nhats(sb_iter) = max(g);

			already_have_eigs = true; % Avoid recomputing eigs
		end%sb_iter

		disp(sb_results);
		for i=1:nresults
			results{i}(pcount,scount) = ops{i}(sb_results{i});
		end
		NHAT(pcount,nsolve) = median(nhats);

		Qtags(pcount) = Qtag;

	end%solvers
end%problems

for i=1:nresults
	eval([result_names{i} '=results{i}']);
end
Q(:,nsolve+1) = Qtags; % gt values

clear benchmark_flag
