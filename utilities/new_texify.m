list_list;
solver_lists{1} = igraph_list;
solver_lists{2} = CV_list;;
problem_lists{1} = six_list;
problem_lists{2} = fb_list;
problem_lists{3} = bio_list;
problem_lists{4} = hsi_list;
nsolver_lists = length(solver_lists);
nproblem_lists = length(problem_lists);

fnum = 1;

for slistnum = 1:nsolver_lists % Solver groups
	solver_list = solver_lists{slistnum};
	nsolvers = length(solver_list);
	for plistnum = 1:nproblem_lists  % Problem groups
		problem_list = problem_lists{plistnum};
		nproblems = length(problem_list);

		if exist('Q'); clear Q; end;
		if exist('Qs'); clear Qs; end;
		if exist('T'); clear T; end;
		if exist('Ts'); clear Ts; end;
		if exist('Te'); clear Te; end;
		if exist('Tes'); clear Tes; end;
		if exist('NMI'); clear NMI; end;
		if exist('NMIs'); clear NMIs; end;
		if exist('PUR'); clear PUR; end;
		if exist('PURs'); clear PURs; end;

		Q_ref_found = false;

		for snum=1:nsolvers; % Solvers
			solver = solver_list{snum};
			for pnum=1:nproblems % Problems
				problem = problem_list{pnum};

				fname = sprintf('%s%s.mat',problem,solver);
				if exist(fname)
					load(fname,'Q','T','Te','NMI','PUR');

					Qs(pnum,snum) = Q(1);
					Ts(pnum,snum) = T;
					Tes(pnum,snum) = Te;
					NMIs(pnum,snum) = NMI;
					PURs(pnum,snum) = PUR;

					Qs(pnum,nsolvers+1) = Q(2); % Get the reference partition value
					Q_ref_found=true;
				else
					disp(['missing file ' fname]);

					Qs(pnum,snum) = 0;
					Ts(pnum,snum) = 0;
					Tes(pnum,snum) = 0;
					NMIs(pnum,snum) = 0;
					PURs(pnum,snum) = 0;

					if ~Q_ref_found
						Qs(pnum,nsolvers+1) = -9999;
					end

				end
			end
		end

		Q = Qs;
		T = Ts;
		Te = Tes;
		NMI = NMIs;
		PUR = PURs;

		cd('benchmarks')
		save(sprintf('test%d.mat',fnum),'Q','T','problem_list','solver_list','Te','NMI','PUR');
		cd('..')
		fnum = fnum+1;
	end
end
